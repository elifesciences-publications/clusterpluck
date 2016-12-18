#!/usr/bin/env Python

import argparse
import sys
import pandas as pd
import numpy as np
import warnings
from clusterpluck.scripts.cluster_dictionary import build_cluster_map
from clusterpluck.scripts.orfs_in_common import generate_index_list
from clusterpluck.scripts.orfs_in_common import pick_a_cluster
from functools import partial
from scoop import futures


# The arg parser
def make_arg_parser():
	parser = argparse.ArgumentParser(description='Report how many ORFs are matched between two clusters with some identity')
	parser.add_argument('-i', '--input', help='Input is the ORF matrix CSV file.', default='-')
	parser.add_argument('-m', '--mpfa', help='The multi-protein fasta file (.mpfa) from which to build the dictionary', required=True)
	parser.add_argument('-b', '--bread', help='Where to find the cluster information in the header for the sequence (default="ref|test,|")', default='ref|,|')
	parser.add_argument('-o', '--output', help='Where to save the output csv; default to screen', required=False, default='-')
	parser.add_argument('-t', '--tanimoto', help='Use the Tanimoto coefficient for the commonality algorithm.', action='store_true', default=False)
	return parser


def generate_chunk_list(in_csv2):
	header = pd.read_csv(in_csv2, header=0, engine='c', index_col=0, nrows=0)
	header = list(header.columns)
	print('Extracted headers from input file...\n')
	return header


def parallel_minicluster(mx, args_list):
	# parse the argument list
	cluster_map = args_list[0]
	c_list = args_list[1]
	cluster = mx.columns[0]
	cluster = '_'.join(cluster.split('|')[3].split('_')[:3])
	j_orfs = len(cluster_map[cluster])  # how many orfs in the full cluster
	c_i = 0
	i = len(c_list)
	mat = np.zeros((i, 1))
	for cluster2 in c_list:
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
		del mx_dubsub
		mat[c_i, 0] = orf_rate
		c_i += 1
	del mx
	print('Finished processing a cluster set!')
	return pd.DataFrame(mat)


def parallel_tanimoto(mx, args_list):
	# parse the argument list
	cluster_map = args_list[0]
	c_list = args_list[1]
	cluster = mx.columns[0]
	cluster = '_'.join(cluster.split('|')[3].split('_')[:3])
	j_orfs = len(cluster_map[cluster])  # how many orfs in the full cluster
	c_i = 0
	i = len(c_list)
	mat = np.zeros((i, 1))
	for cluster2 in c_list:
		i_orfs = len(cluster_map[cluster2])
		# subsets the smaller matrix by rows belonging to one cluster
		mx = pd.DataFrame(mx)
		mx_dubsub = mx.filter(like=cluster2, axis=0)
		mx_dubsub = mx_dubsub.dropna(axis=(1, 0), how='all')
		i_mx = mx_dubsub.shape[0]
		j_mx = mx_dubsub.shape[1]
		nc = min([i_mx, j_mx])
		with warnings.catch_warnings():
			warnings.simplefilter('ignore', category=RuntimeWarning)  # np doesn't like taking mean of empty slices
			cc_mean = np.nanmean(mx_dubsub.values, dtype='float64')
		if cc_mean > 0:
			# calculates the Tanimoto coefficient for the orf-orf comparison
			tanimoto = nc / (i_orfs + j_orfs - nc)
		else:
			tanimoto = cc_mean
		# saves this ratio into the pre-existing array at the (cluster, cluster2) location
		del mx_dubsub
		mat[c_i, 0] = tanimoto
		c_i += 1
	del mx
	print('Finished processing a cluster set!')
	return pd.DataFrame(mat)


def main():
	parser = make_arg_parser()
	args = parser.parse_args()
	# Parse command line
	tanimoto = args.tanimoto
	with open(args.mpfa, 'r') as inf:
		# Generates dictionary with each unique 'refseq_cluster' as keys, ORFs as values
		cluster_map = build_cluster_map(inf, bread=args.bread)
	with open(args.input, 'r') as inf2:
		inkey = generate_index_list(inf2)
		print('\nOk, processing input file...\n')
	with open(args.input, 'r') as in_csv2:
		headers = generate_chunk_list(in_csv2)
	c_list = list(cluster_map.keys())
	grabbed_clusters = []
	data_to_pool = []
	# print(c_list)
	for cluster in c_list:
		grab = pick_a_cluster(headers, cluster)  # uses the name of the cluster to get a list of all orfs for a particular unique cluster
		# print(grab)
		if not grab:
			pass
		else:
			# print(grab)
			grabbed_clusters.extend([cluster])
			with open(args.input, 'r') as inf3:
				mx = pd.read_csv(inf3, sep=',', header=0, usecols=grab, engine='c')  # loads in only the columns from the grab list, i.e. all cols for a unique cluster
			mx.index = inkey  # reindexes the df with the orf labels after importing specific columns with usecols
			data_to_pool.append(mx)
	# print(grabbed_clusters)
	# print(len(data_to_pool))
	args_list = [cluster_map, c_list]  # organizes all the arguments that the parallelized function needs into a list
	if __name__ == '__main__':
		print('\nSending data to Workers... work, Workers, work!\n')
		if tanimoto:
			results = list(futures.map(partial(parallel_tanimoto, args_list=args_list), data_to_pool))
		else:
			results = list(futures.map(partial(parallel_minicluster, args_list=args_list), data_to_pool))
		# bigmat = pd.concat(results, axis=0)  # stack all the results into a single column in a dataframe
		# print(bigmat.shape[0])
		# bigmat.index = c_list  # now the index is just the clusters, not the orfs
		# print(bigmat)
	print('File processing complete; writing output file...\n')
	with open(args.output, 'w') if args.output != '-' else sys.stdout as outf:
		del data_to_pool
		outdf = pd.concat(results, axis=1)
		outdf.columns = grabbed_clusters  # names the columns (and index, next line) according to clusters in the order they were processed
		outdf.index = c_list
		outdf.sort_index(axis=0, inplace=True)
		outdf.sort_index(axis=1, inplace=True)
		outdf = outdf.round(decimals=3)
		outdf.to_csv(outf)


if __name__ == '__main__':
	main()
