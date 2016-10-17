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


def cluster_completeness(cluster_map, in_csv):
	mx = pd.read_csv(in_csv, sep=',', header=0, index_col=0)
	# list of all the clusters
	cluster_coverage = defaultdict(dict)
	for cluster in cluster_map:
		# how many orfs in the full cluster
		j_orfs = len(cluster_map[cluster])
		# subsets the matrix by columns belonging to one cluster
		mx_csub = mx.filter(like=cluster)
		for cluster2 in cluster_map:
			i_orfs = len(cluster_map[cluster2])
			# subsets the smaller matrix by rows belonging to one cluster
			mx_dubsub = mx_csub.filter(like=cluster2, axis=0)
			mx_dubsub.dropna(axis=(1, 0), how='all', inplace=True)
			i_mx = mx_dubsub.shape[0]
			j_mx = mx_dubsub.shape[1]
			cc_mean = np.nanmean(mx_dubsub.values, dtype='float64')
			if cc_mean > 0:
				# calculates the fraction of orf coverage by the match
				orf_rate = (j_mx + i_mx) / (j_orfs + i_orfs)  # one way of doing it
				# orf_rate = ((j_mx / j_orfs) + (i_mx / i_orfs)) / 2  # alternative metric
			else:
				orf_rate = cc_mean
			# saves this ratio in a dictionary
			cluster_coverage[cluster][cluster2] = orf_rate
	coverage_matrix = pd.DataFrame.from_dict(cluster_coverage, orient='columns', dtype=float)
	coverage_matrix.sort_index(axis=0, inplace=True)
	coverage_matrix.sort_index(axis=1, inplace=True)
	# print(coverage_matrix)
	return coverage_matrix


def main():
	parser = make_arg_parser()
	args = parser.parse_args()
	# Parse command line
	with open(args.mpfa, 'r') if args.mpfa != '-' else sys.stdin as inf:
		# Generates dictionary with each unique 'refseq_cluster' as keys, ORFs as values
		cluster_map = build_cluster_map(inf, bread=args.bread)
		with open(args.input, 'r') as in_csv:
			with open(args.output, 'w') if args.output != '-' else sys.stdout as outf:
				coverage_matrix = cluster_completeness(cluster_map, in_csv)
				coverage_matrix = coverage_matrix.round(decimals=3)
				coverage_matrix.to_csv(outf)

if __name__ == '__main__':
	with warnings.catch_warnings():
		warnings.filterwarnings('ignore', message=r'Mean of empty slice', category=RuntimeWarning)
		main()
