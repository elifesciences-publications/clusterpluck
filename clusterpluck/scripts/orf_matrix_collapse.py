#!/usr/bin/env Python

import argparse
import sys
import numpy as np
import pandas as pd
import warnings

from collections import defaultdict
from clusterpluck.scripts.cluster_dictionary import build_cluster_map


# The arg parser
def make_arg_parser():
	parser = argparse.ArgumentParser(description='collapse the ORF x ORF matrix into a scored cluster x cluster matrix')
	parser.add_argument('-i', '--input', help='Input is the ORF matrix CSV file.', default='-')
	parser.add_argument('-m', '--mpfa', help='The multi-protein fasta file (.mpfa) from which to build the dictionary')
	parser.add_argument('-b', '--bread', help='Where to find the cluster information in the header for the sequence (default="ref|,|")', default='ref|,|')
	parser.add_argument('-o', '--output', help='Where to save the output csv; default to screen', required=False, default='-')
	return parser


def cluster_by_cluster(cluster_map, in_csv):
	mx = pd.read_csv(in_csv, sep=',', header=0, index_col=0)
	c_list = list(cluster_map.keys())
	ct = len(c_list)
	mat = np.zeros((ct, ct))
	j = 0
	for cluster in cluster_map:
		# subsets the matrix by columns belonging to one cluster
		mx_csub = mx.filter(like=cluster)
		i = 0
		for cluster2 in cluster_map:
			# subsets the smaller matrix by rows belonging to one cluster
			mx_dubsub = mx_csub.filter(like=cluster2, axis=0)
			# finds the mean of the cells in the cluster x cluster matrix
			with warnings.catch_warnings():
				warnings.simplefilter('ignore', category=RuntimeWarning)
				cc_mean = np.nanmean(mx_dubsub.values, dtype='float64')
			# saves this average in a dictionary
			mat[i, j] = cc_mean
			i += 1
		j += 1
	outdf = pd.DataFrame(mat, dtype=float)
	outdf.columns = c_list
	outdf.index = c_list
	outdf.sort_index(axis=0)
	outdf.sort_index(axis=1)
	# print(score_mean)
	return outdf
	# Check if a matrix is symmetric
	# arr = df.values
	# print((arr.transpose() == -arr).all())


def main():
	parser = make_arg_parser()
	args = parser.parse_args()
	# Parse command line
	with open(args.mpfa, 'r') if args.mpfa != '-' else sys.stdin as inf:
		# Generates dictionary with each unique 'refseq_cluster' as keys, ORFs as values
		cluster_map = build_cluster_map(inf, bread=args.bread)
		with open(args.input, 'r') as in_csv:
			with open(args.output, 'w') if args.output != '-' else sys.stdout as outf:
				outdf = cluster_by_cluster(cluster_map, in_csv)
				outdf = outdf.round(decimals=2)
				outdf.to_csv(outf)

if __name__ == '__main__':
	with warnings.catch_warnings():
		warnings.filterwarnings('ignore')
		main()
