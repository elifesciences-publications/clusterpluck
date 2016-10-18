#!/usr/bin/env Python

import argparse
import sys
import pandas as pd
import numpy as np
import warnings
from collections import defaultdict
from clusterpluck.scripts.cluster_dictionary import build_cluster_map


# The arg parser
def make_arg_parser():
	parser = argparse.ArgumentParser(description='Report how many ORFs are matched between two clusters with some identity')
	parser.add_argument('-i', '--input', help='Input is the ORF matrix CSV file.', default='-')
	parser.add_argument('-m', '--mpfa', help='The multi-protein fasta file (.mpfa) from which to build the dictionary')
	parser.add_argument('-b', '--bread', help='Where to find the cluster information in the header for the sequence (default="ref|,|")', default='ref|,|')
	parser.add_argument('-o', '--output', help='Where to save the output csv; default to screen', required=False, default='-')
	return parser


def cluster_completeness(intype, cluster_map, inf2):
	if intype == 'csv':
		mx = pd.read_csv(inf2, sep=',', header=0, index_col=0)
	elif intype == 'h5':
		mx = pd.read_hdf(inf2, 'table')
	# list of all the clusters
	# cluster_coverage = defaultdict(dict)
	c_list = list(cluster_map.keys())
	ct = len(c_list)
	mat = np.zeros((ct, ct))
	j = 0
	for cluster in c_list:
		# print(cluster)
		# how many orfs in the full cluster
		j_orfs = len(cluster_map[cluster])
		# subsets the matrix by columns belonging to one cluster
		mx_csub = mx.filter(like=cluster)
		i = 0
		for cluster2 in c_list:
			# print(cluster2)
			# print(i)
			# print(j)
			i_orfs = len(cluster_map[cluster2])
			# subsets the smaller matrix by rows belonging to one cluster
			mx_dubsub = mx_csub.filter(like=cluster2, axis=0)
			mx_dubsub.dropna(axis=(1, 0), how='all', inplace=True)
			i_mx = mx_dubsub.shape[0]
			j_mx = mx_dubsub.shape[1]
			with warnings.catch_warnings():
				warnings.simplefilter('ignore', category=RuntimeWarning)
				cc_mean = np.nanmean(mx_dubsub.values, dtype='float64')
			if cc_mean > 0:
				# calculates the fraction of orf coverage by the match
				orf_rate = (j_mx + i_mx) / (j_orfs + i_orfs)  # one way of doing it
				# orf_rate = ((j_mx / j_orfs) + (i_mx / i_orfs)) / 2  # alternative metric
			else:
				orf_rate = cc_mean
			# saves this ratio in a dictionary
			mat[i, j] = orf_rate
			i += 1
		j += 1
	outdf = pd.DataFrame(mat, dtype=float)
	outdf.columns = c_list
	outdf.index = c_list
	outdf.sort_index(axis=0, inplace=True)
	outdf.sort_index(axis=1, inplace=True)
	# print(outdf)
	return outdf


def main():
	parser = make_arg_parser()
	args = parser.parse_args()
	# Parse command line
	with open(args.mpfa, 'r') if args.mpfa != '-' else sys.stdin as inf:
		# Generates dictionary with each unique 'refseq_cluster' as keys, ORFs as values
		cluster_map = build_cluster_map(inf, bread=args.bread)
		intype = str(args.input).split('.')[-1]
		with open(args.input, 'r') as inf2:
			with open(args.output, 'w') if args.output != '-' else sys.stdout as outf:
				outdf = cluster_completeness(intype, cluster_map, inf2)
				outdf = outdf.round(decimals=3)
				outdf.to_csv(outf)

if __name__ == '__main__':
	main()
