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

# usage = python -m scoop -vv -n 480 mpi_collapse -i [input.csv] -m [.mpfa key] -o [output.csv]

# The arg parser
def make_arg_parser():
	parser = argparse.ArgumentParser(description='Collapse ORF matrix into a scored cluster matrix. Run with "python -m scoop -vv -n 480 mpi_parallel_collapse [args]"')
	parser.add_argument('-i', '--input', help='Input is the ORF matrix CSV file.', default='-')
	parser.add_argument('-m', '--mpfa', help='The multi-protein fasta file (.mpfa) from which to build the dictionary')
	parser.add_argument('-b', '--bread', help='Where to find the cluster information in the header for the sequence (default="ref|,|")', default='ref|,|')
	parser.add_argument('-o', '--output', help='Where to save the output csv; default to screen', required=False, default='-')
	return parser


def generate_chunk_list(in_csv2):
	header = pd.read_csv(in_csv2, header=0, engine='c', index_col=0, nrows=0)
	header = list(header.columns)
	print('Extracted headers from input file...\n')
	return header


def parallel_clustermean(mx, c_list):
	i = len(c_list)
	mat = np.zeros((i, 1))
	c_i = 0
	for cluster2 in c_list:
		mx_dubsub = mx.filter(like=cluster2, axis=0)  # subsets the smaller matrix by rows belonging to one cluster
		# finds the mean of the cells in the cluster x cluster2 matrix
		with warnings.catch_warnings():
			warnings.simplefilter('ignore', category=RuntimeWarning)  # np doesn't like taking mean of empty slices
			cc_mean = np.nanmean(mx_dubsub.values, dtype='float64')
		mat[c_i, 0] = cc_mean  # saves this mean into the pre-existing array at the right location
		c_i += 1
	del mx
	dfmeans = pd.DataFrame(mat)
	dfmeans = dfmeans.round(decimals=2)
	dfmeans.index = c_list
	return dfmeans


def main():
	parser = make_arg_parser()
	args = parser.parse_args()
	# Parse command line
	with open(args.mpfa, 'r') as inf:
		# Generates dictionary with each unique 'refseq_cluster' as keys, ORFs as values
		cluster_map = build_cluster_map(inf, bread=args.bread)
	with open(args.input, 'r') as in_csv:
		print('\nOk, processing input file in pieces...\n')
		inkey = generate_index_list(in_csv)
		# print(len(inkey))
	with open(args.input, 'r') as in_csv2:
		headers = generate_chunk_list(in_csv2)
		# print(len(headers))
		c_list = list(cluster_map.keys())
		# ct = len(c_list)
		# print('Found %d clusters...' % ct)
		data_to_pool = []
		grabbed_clusters = []
	for cluster in c_list:
		grab = pick_a_cluster(headers, cluster)  # uses the name of the cluster to get a list of all orfs for a particular unique cluster
		if not grab:
			pass
		else:
			# print(grab)
			grabbed_clusters.extend([cluster])
			with open(args.input, 'r') as inf3:
				mx = pd.read_csv(inf3, sep=',', header=0, usecols=grab, engine='c')  # loads in only the columns from the grab list, i.e. all cols for a unique cluster
			mx.index = inkey  # reindexes the df with the orf labels after importing specific columns with usecols
			data_to_pool.append(mx)  # create the list of dfs to map over for multiprocessing
	if __name__ == '__main__':
		print('\nSending data to Workers... work, Workers, work!')
		results = list(futures.map(partial(parallel_clustermean, c_list=c_list), data_to_pool))
		print('\nFile processing complete; writing output file...\n')
		del data_to_pool
	with open(args.output, 'w') if args.output != '-' else sys.stdout as outf:
		outdf = pd.concat(results, axis=1)
		outdf.columns = grabbed_clusters  # names the columns (and index, next line) according to clusters in the order they were processed
		# outdf.index = c_list
		outdf.sort_index(axis=0, inplace=True)  # ensure that the clusters are in order on cols and rows
		outdf.sort_index(axis=1, inplace=True)
		outdf.to_csv(outf)


if __name__ == '__main__':
	main()
