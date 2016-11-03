#!/usr/bin/env Python

import argparse
import sys
import pandas as pd
import numpy as np
import warnings
import os
from clusterpluck.scripts.cluster_dictionary import build_cluster_map
from functools import partial
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
	return parser


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


def big_cluster_completeness(inf3, grab, inkey, cluster, cluster_map, cpus, c_list, j):
	mx = pd.read_csv(inf3, sep=',', header=0, usecols=grab, engine='c')
	mx.index = inkey
	# how many orfs in the full cluster
	j_orfs = len(cluster_map[cluster])
	args_list = [mx, j_orfs, cluster_map]
	with Pool(processes=cpus) as pool:
		results = pool.starmap(parallel_minicluster, zip(c_list, repeat(args_list)))
		pool.close()
		pool.join()
		bigmat = pd.concat(results, axis=0)
		# print(bigmat.shape[0])
		bigmat.index = c_list
	# DEBUG - will print the progress every 100 clusters (across the slower dimension).
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
	if cpus > num_cpus:
		print('\nError: Number of requested processors exceeds hardware available!')
		print('Maximum processors available is %s.\n' % num_cpus)
		sys.exit()
	with open(args.mpfa, 'r') as inf:
		# Generates dictionary with each unique 'refseq_cluster' as keys, ORFs as values
		cluster_map = build_cluster_map(inf, bread=args.bread)
		intype = str(args.input).split('.')[-1]
		with open(args.input, 'r') as inf2:
			with open(args.output, 'w') if args.output != '-' else sys.stdout as outf:
				print('\nOk, processing input file in pieces...\n')
				inkey = generate_index_list(inf2)
				c_list = list(cluster_map.keys())
				ct = len(c_list)
				print('Found %d clusters...' % ct)
				results_list = []
				j = 0
				for cluster in c_list:
					grab = pick_a_cluster(inkey, cluster)
					# print(grab)
					with open(args.input, 'r') as inf3:
						bigmat = big_cluster_completeness(inf3, grab, inkey, cluster, cluster_map, cpus, c_list, j)
					# print(bigmat)
					results_list.append(bigmat)
					j += 1
				print('File processing complete; writing output file...\n')
				outdf = pd.concat(results_list, axis=1)
				outdf.columns = c_list  # names the columns (and index, next line) according to clusters in the order they were processed
				outdf.index = c_list
				outdf.sort_index(axis=0, inplace=True)
				outdf.sort_index(axis=1, inplace=True)
				outdf = outdf.round(decimals=3)
				outdf.to_csv(outf)


if __name__ == '__main__':
	main()
