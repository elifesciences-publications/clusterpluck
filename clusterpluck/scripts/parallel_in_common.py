#!/usr/bin/env Python

import argparse
import sys
import pandas as pd
import numpy as np
import warnings
from clusterpluck.scripts.cluster_dictionary import build_cluster_map
from itertools import repeat
from multiprocessing import Pool
from multiprocessing import cpu_count


# The arg parser
def make_arg_parser():
	parser = argparse.ArgumentParser(description='Report how many ORFs are matched between two clusters with some identity')
	parser.add_argument('-i', '--input', help='Input is the ORF matrix CSV file.', default='-')
	parser.add_argument('-m', '--mpfa', help='The multi-protein fasta file (.mpfa) from which to build the dictionary', required=True)
	parser.add_argument('-b', '--bread', help='Where to find the cluster information in the header for the sequence (default="ref|test,|")', default='ref|,|')
	parser.add_argument('-o', '--output', help='Where to save the output csv; default to screen', required=False, default='-')
	parser.add_argument('-c', '--cpus', help='How many processors to use (integer); default is 4', required=False, default=4, type=int)
	parser.add_argument('-t', '--tanimoto', help='Use the Tanimoto coefficient for the commonality algorithm.', action='store_true', default=False)
	return parser


def generate_chunk_list(in_csv2):
	header = pd.read_csv(in_csv2, header=0, engine='c', index_col=0, nrows=0)
	header = list(header.columns)
	print('Extracted headers from input file...\n')
	return header


def generate_index_list(inf2):
	rowdf = pd.read_csv(inf2, header=0, engine='c', usecols=[0])
	inkey = rowdf.iloc[0:, 0].tolist()
	del rowdf
	print('Extracted cluster names from input file...\n')
	return inkey


def pick_a_cluster(inkey, cluster):
	grab = [n for n in inkey if cluster in n]
	return grab


def parallel_minicluster(cluster2, args_list):
	# parse the argument list
	mx = args_list[0]
	j_orfs = args_list[1]
	cluster_map = args_list[2]
	i_orfs = len(cluster_map[cluster2])
	# subsets the smaller matrix by rows belonging to one cluster
	mx_dubsub = mx.filter(like=cluster2, axis=0)
	mx_dubsub = mx_dubsub.dropna(axis=(1, 0), how='all')
	i_mx = mx_dubsub.shape[0]
	j_mx = mx_dubsub.shape[1]
	with warnings.catch_warnings():
		warnings.simplefilter('ignore', category=RuntimeWarning)  # np doesn't like taking mean of empty slices
		cc_mean = np.nanmean(mx_dubsub.values, dtype='float64')
	if cc_mean > 0:
		# calculates the fraction of orf coverage by the match
		orf_rate = (j_mx + i_mx) / (j_orfs + i_orfs)  # one way of doing it
	# orf_rate = ((j_mx / j_orfs) + (i_mx / i_orfs)) / 2  # alternative metric
	else:
		orf_rate = cc_mean
	# saves this ratio into the pre-existing array at the (cluster, cluster2) location
	mat = np.zeros((1, 1))
	mat[0, 0] = orf_rate
	return pd.DataFrame(mat)


def parallel_tanimoto(cluster2, args_list):
	# parse the argument list
	mx = args_list[0]
	j_orfs = args_list[1]
	cluster_map = args_list[2]
	i_orfs = len(cluster_map[cluster2])
	# subsets the smaller matrix by rows belonging to one cluster
	mx_dubsub = mx.filter(like=cluster2, axis=0)
	mx_dubsub = mx_dubsub.dropna(axis=(1, 0), how='all')
	i_mx = mx_dubsub.shape[0]
	j_mx = mx_dubsub.shape[1]
	nc = min([i_mx, j_mx])
	with warnings.catch_warnings():
		warnings.simplefilter('ignore', category=RuntimeWarning)  # np doesn't like taking mean of empty slices
		cc_mean = np.nanmean(mx_dubsub.values, dtype='float64')
	if cc_mean > 0:
		# calculates the fraction of orf coverage by the match
		tanimoto = nc / (i_orfs + j_orfs - nc)
	else:
		tanimoto = cc_mean
	# saves this ratio into the pre-existing array at the (cluster, cluster2) location
	mat = np.zeros((1, 1))
	mat[0, 0] = tanimoto
	return pd.DataFrame(mat)


def big_cluster_completeness(inf3, grab, inkey, cluster, cluster_map, cpus, tanimoto, c_list, j):
	mx = pd.read_csv(inf3, sep=',', header=0, usecols=grab, engine='c')  # loads in only the columns from the grab list, i.e. all cols for a unique cluster
	mx.index = inkey  # reindexes the df with the orf labels after importing specific columns with usecols
	# how many orfs in the full cluster
	j_orfs = len(cluster_map[cluster])
	args_list = [mx, j_orfs, cluster_map]  # organizes all the arguments that the parallelized function needs into a list
	with Pool(processes=cpus) as pool:
		if not tanimoto:
			results = pool.starmap(parallel_minicluster, zip(c_list, repeat(args_list)))
		else:
			results = pool.starmap(parallel_tanimoto, zip(c_list, repeat(args_list)))
		pool.close()
		pool.join()
		bigmat = pd.concat(results, axis=0)  # stack all the results into a single column in a dataframe
		# print(bigmat.shape[0])
		bigmat.index = c_list  # now the index is just the clusters, not the orfs
	# DEBUG - will print the progress every 50 clusters (across columns--the slower dimension).
	if j % 50:
		pass
	elif j == 0:
		print('Processed first cluster... moving on!')
	else:
		print('Processed %d clusters' % j)
	del mx
	return bigmat


def main():
	parser = make_arg_parser()
	args = parser.parse_args()
	# Parse command line
	cpus = args.cpus
	num_cpus = cpu_count()
	tanimoto = args.tanimoto
	if cpus > num_cpus:
		print('\nError: Number of requested processors exceeds hardware available!')
		print('Maximum processors available is %s.\n' % num_cpus)
		sys.exit()
	with open(args.mpfa, 'r') as inf:
		# Generates dictionary with each unique 'refseq_cluster' as keys, ORFs as values
		cluster_map = build_cluster_map(inf, bread=args.bread)
		# intype = str(args.input).split('.')[-1]
	with open(args.input, 'r') as inf2:
		inkey = generate_index_list(inf2)
	with open(args.input, 'r') as in_csv2:
		headers = generate_chunk_list(in_csv2)
	c_list = list(cluster_map.keys())
	ct = len(c_list)
	print('Found %d clusters...' % ct)
	results_list = []
	grabbed_clusters = []
	j = 0
	for cluster in c_list:
		grab = pick_a_cluster(headers, cluster)  # uses the name of the cluster to get a list of all orfs for a particular unique cluster
		# print(grab)
		if not grab:
			pass
		else:
			grabbed_clusters.extend([cluster])
			with open(args.input, 'r') as inf3:
				bigmat = big_cluster_completeness(inf3, grab, inkey, cluster, cluster_map, cpus, tanimoto, c_list, j)
			results_list.append(bigmat)  # returns a list of dataframes, one for each cluster column
			j += 1
	print('File processing complete; writing output file...\n')
	with open(args.output, 'w') if args.output != '-' else sys.stdout as outf:
		outdf = pd.concat(results_list, axis=1)
		outdf.columns = grabbed_clusters  # names the columns (and index, next line) according to clusters in the order they were processed
		outdf.index = c_list
		outdf.sort_index(axis=0, inplace=True)
		outdf.sort_index(axis=1, inplace=True)
		outdf = outdf.round(decimals=3)
		outdf.to_csv(outf)


if __name__ == '__main__':
	main()
