#!/usr/bin/env Python

import argparse
import sys
import numpy as np
import pandas as pd
import warnings
from clusterpluck.scripts.cluster_dictionary import build_cluster_map
from clusterpluck.scripts.orfs_in_common import generate_index_list
from clusterpluck.scripts.orfs_in_common import pick_a_cluster
from functools import partial
from scoop import futures


# The arg parser
def make_arg_parser():
	parser = argparse.ArgumentParser(description='Collapse ORF matrix into a scored cluster matrix. Run with "python -m scoop -vv -n 480 mpi_parallel_collapse [args]"')
	parser.add_argument('-i', '--input', help='Input is the ORF matrix CSV file.', default='-')
	parser.add_argument('-m', '--mpfa', help='The multi-protein fasta file (.mpfa) from which to build the dictionary')
	parser.add_argument('-b', '--bread', help='Where to find the cluster information in the header for the sequence (default="ref|,|")', default='ref|,|')
	parser.add_argument('-o', '--output', help='Where to save the output csv; default to screen', required=False, default='-')
	return parser


def parallel_clustermean(cluster2, mx):
	# subsets the smaller matrix by rows belonging to one cluster
	mx_dubsub = mx.filter(like=cluster2, axis=0)
	# finds the mean of the cells in the cluster x cluster2 matrix
	with warnings.catch_warnings():
		warnings.simplefilter('ignore', category=RuntimeWarning)  # np doesn't like taking mean of empty slices
		cc_mean = np.nanmean(mx_dubsub.values, dtype='float64')
	# saves this ratio into the pre-existing array at the (cluster, cluster2) location
	mat = np.zeros((1, 1))
	mat[0, 0] = cc_mean
	return pd.DataFrame(mat)


def big_cluster_v_cluster(inf3, grab, inkey, c_list, j):
	mx = pd.read_csv(inf3, sep=',', header=0, usecols=grab, engine='c')  # loads in only the columns from the grab list, i.e. all cols for a unique cluster
	mx.index = inkey  # reindexes the df with the orf labels after importing specific columns with usecols
	# how many orfs in the full cluster
	# args_list = [mx]  # organizes all the arguments that the parallelized function needs into a list
	if __name__ == '__main__':
		results = list(futures.map(partial(parallel_clustermean, mx=mx), c_list))
		bigmat = pd.concat(results, axis=0)  # stack all the results into a single column in a dataframe
		# print(bigmat.shape[0])
		bigmat.index = c_list  # now the index is just the clusters, not the orfs
	# DEBUG - will print the progress every 50 clusters (across the slower dimension).
	if j % 50:
		pass
	elif j == 0:
		print('Processed first cluster... moving on!')
	elif j == 1:
		print('Worry not, second cluster has been processed...')
	else:
		print('Processed %d clusters' % j)
	del mx
	return bigmat


def main():
	parser = make_arg_parser()
	args = parser.parse_args()
	# Parse command line
	with open(args.mpfa, 'r') as inf:
		# Generates dictionary with each unique 'refseq_cluster' as keys, ORFs as values
		cluster_map = build_cluster_map(inf, bread=args.bread)
	with open(args.input, 'r') as in_csv:
		with open(args.output, 'w') if args.output != '-' else sys.stdout as outf:
			print('\nOk, processing input file in pieces...\n')
			inkey = generate_index_list(in_csv)
			c_list = list(cluster_map.keys())
			ct = len(c_list)
			print('Found %d clusters...' % ct)
			results_list = []
			j = 0
			for cluster in c_list:
				grab = pick_a_cluster(inkey, cluster)  # uses the name of the cluster to get a list of all orfs for a particular unique cluster
				# print(grab)
				with open(args.input, 'r') as inf3:
					bigmat = big_cluster_v_cluster(inf3, grab, inkey, c_list, j)
				# print(bigmat)
				results_list.append(bigmat)  # returns a list of dataframes, one for each cluster column
				j += 1
			print('File processing complete; writing output file...\n')
			outdf = pd.concat(results_list, axis=1)
			outdf.columns = c_list  # names the columns (and index, next line) according to clusters in the order they were processed
			outdf.index = c_list
			outdf.sort_index(axis=0, inplace=True)
			outdf.sort_index(axis=1, inplace=True)
			outdf = outdf.round(decimals=2)
			outdf.to_csv(outf)


if __name__ == '__main__':
	main()
