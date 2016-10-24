#!/usr/bin/env Python

import argparse
import sys
import numpy as np
import pandas as pd
import warnings
from clusterpluck.scripts.cluster_dictionary import build_cluster_map
from clusterpluck.scripts.orfs_in_common import generate_index_list
from clusterpluck.scripts.orfs_in_common import pick_a_cluster


# The arg parser
def make_arg_parser():
	parser = argparse.ArgumentParser(description='collapse the ORF x ORF matrix into a scored cluster x cluster matrix')
	parser.add_argument('-i', '--input', help='Input is the ORF matrix CSV file.', default='-')
	parser.add_argument('-m', '--mpfa', help='The multi-protein fasta file (.mpfa) from which to build the dictionary')
	parser.add_argument('-b', '--bread', help='Where to find the cluster information in the header for the sequence (default="ref|,|")', default='ref|,|')
	parser.add_argument('-o', '--output', help='Where to save the output csv; default to screen', required=False, default='-')
	parser.add_argument('-p', '--pieces', help='Read the input CSV in pieces, useful for very large (>5GB) files', action='store_true', required=False, default=False)
	return parser


def cluster_by_cluster(cluster_map, in_csv):
	mx = pd.read_csv(in_csv, sep=',', header=0, index_col=0)
	c_list = list(cluster_map.keys())  # list of all clusters
	ct = len(c_list)
	mat = np.zeros((ct, ct))  # initializes an array to fit results from all clusters
	j = 0
	for cluster in cluster_map:
		# subsets the matrix by columns belonging to one cluster
		mx_csub = mx.filter(like=cluster)
		i = 0
		for cluster2 in cluster_map:
			# subsets the smaller matrix by rows belonging to one cluster
			mx_dubsub = mx_csub.filter(like=cluster2, axis=0)
			# finds the mean of the cells in the cluster x cluster2 matrix
			with warnings.catch_warnings():
				warnings.simplefilter('ignore', category=RuntimeWarning)  # np doesn't like taking mean of empty slices
				cc_mean = np.nanmean(mx_dubsub.values, dtype='float64')
			# saves this ratio into the pre-existing array at the (cluster, cluster2) location
			mat[i, j] = cc_mean
			i += 1
		j += 1
	outdf = pd.DataFrame(mat, dtype=float)
	outdf.columns = c_list  # names the columns (and index, next line) according to clusters in the order they were processed
	outdf.index = c_list
	outdf.sort_index(axis=0, inplace=True)
	outdf.sort_index(axis=1, inplace=True)
	# print(score_mean)
	return outdf
	# Check if a matrix is symmetric
	# arr = df.values
	# print((arr.transpose() == -arr).all())


def big_cluster_by_cluster(grab, inkey, c_list, inf3, mat, j):
	mx = pd.read_csv(inf3, sep=',', header=0, usecols=grab, engine='c')
	mx.index = inkey
	# how many orfs in the full cluster
	i = 0
	for cluster2 in c_list:
		# subsets the smaller matrix by rows belonging to one cluster
		mx_dubsub = mx.filter(like=cluster2, axis=0)
		# finds the mean of the cells in the cluster x cluster2 matrix
		with warnings.catch_warnings():
			warnings.simplefilter('ignore', category=RuntimeWarning)  # np doesn't like taking mean of empty slices
			cc_mean = np.nanmean(mx_dubsub.values, dtype='float64')
		# saves this ratio into the pre-existing array at the (cluster, cluster2) location
		mat[i, j] = cc_mean
		i += 1
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
	return mat


def main():
	parser = make_arg_parser()
	args = parser.parse_args()
	# Parse command line
	with open(args.mpfa, 'r') if args.mpfa != '-' else sys.stdin as inf:
		# Generates dictionary with each unique 'refseq_cluster' as keys, ORFs as values
		cluster_map = build_cluster_map(inf, bread=args.bread)
		with open(args.input, 'r') as in_csv:
			with open(args.output, 'w') if args.output != '-' else sys.stdout as outf:
				if not args.pieces:
					outdf = cluster_by_cluster(cluster_map, in_csv)
					outdf = outdf.round(decimals=2)
					outdf.to_csv(outf)
				else:
					print('\nOk, processing input file in pieces...\n')
					inkey = generate_index_list(in_csv)
					c_list = list(cluster_map.keys())
					ct = len(c_list)
					print('Found %d clusters...' % ct)
					mat = np.zeros((ct, ct))  # initializes an array of the dimensions necessary to fit all cluster results
					j = 0
					for cluster in c_list:
						grab = pick_a_cluster(inkey, cluster)
						# print(grab)
						with open(args.input, 'r') as inf3:
							mat = big_cluster_by_cluster(grab, inkey, c_list, inf3, mat, j)
						# print(mat)
						j += 1
					print('File processing complete; writing output file...\n')
					outdf = pd.DataFrame(mat, dtype=float)
					outdf.columns = c_list  # names the columns (and index, next line) according to clusters in the order they were processed
					outdf.index = c_list
					outdf.sort_index(axis=0, inplace=True)
					outdf.sort_index(axis=1, inplace=True)
					outdf = outdf.round(decimals=2)
					outdf.to_csv(outf)


if __name__ == '__main__':
	main()
