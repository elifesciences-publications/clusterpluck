#!/usr/bin/env Python

import argparse
import sys
import pandas as pd
import numpy as np
import warnings
import os
from clusterpluck.scripts.cluster_dictionary import build_cluster_map


# The arg parser
def make_arg_parser():
	parser = argparse.ArgumentParser(description='Report how many ORFs are matched between two clusters with some identity')
	parser.add_argument('-i', '--input', help='Input is the ORF matrix CSV file.', default='-')
	parser.add_argument('-m', '--mpfa', help='The multi-protein fasta file (.mpfa) from which to build the dictionary')
	parser.add_argument('-b', '--bread', help='Where to find the cluster information in the header for the sequence (default="ref|,|")', default='ref|,|')
	parser.add_argument('-o', '--output', help='Where to save the output csv; default to screen', required=False, default='-')
	parser.add_argument('-p', '--pieces', help='Read the input CSV in pieces, useful for very large (>5GB) files', action='store_true', required=False, default=False)
	return parser


def cluster_completeness(intype, cluster_map, inf2):
	if intype == 'csv':
		mx = pd.read_csv(inf2, sep=',', header=0, index_col=0, engine='c')
	elif intype == 'h5':
		mx = pd.read_hdf(inf2, 'table')
	# list of all the clusters
	c_list = list(cluster_map.keys())
	ct = len(c_list)
	mat = np.zeros((ct, ct))  # initializes an array of the dimensions necessary to fit all cluster results
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
			mat[i, j] = orf_rate
			i += 1
		j += 1
		# DEBUG - will print the progress every 100 clusters (across the slower dimension) every 500th orf.
		if j % 100:
			pass
		else:
			print(j)
	outdf = pd.DataFrame(mat, dtype=float)
	outdf.columns = c_list  # names the columns (and index, next line) according to clusters in the order they were processed
	outdf.index = c_list
	outdf.sort_index(axis=0, inplace=True)
	outdf.sort_index(axis=1, inplace=True)
	# print(outdf)
	return outdf


def generate_index_list(inf2):
	rowdf = pd.read_csv(inf2, header=0, engine='c', usecols=[0])
	inkey = rowdf.iloc[0:, 0].tolist()
	del rowdf
	print('Extracted cluster names from input file...\n')
	return inkey


def pick_a_cluster(inkey, cluster):
	grab = [n for n in inkey if cluster in n]
	return grab


def big_cluster_completeness(grab, inkey, cluster, cluster_map, c_list, inf3, mat, j):
	mx = pd.read_csv(inf3, sep=',', header=0, usecols=grab, engine='c')
	mx.index = inkey
	# how many orfs in the full cluster
	j_orfs = len(cluster_map[cluster])
	i = 0
	for cluster2 in c_list:
		# print(cluster2)
		# print(i)
		# print(j)
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
		mat[i, j] = orf_rate
		i += 1
	# DEBUG - will print the progress every 100 clusters (across the slower dimension).
	if j % 100:
		pass
	elif j == 0:
		print('Processed first cluster... moving on!')
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
		intype = str(args.input).split('.')[-1]
		insize = os.stat(args.input).st_size
		with open(args.input, 'r') as inf2:
			with open(args.output, 'w') if args.output != '-' else sys.stdout as outf:
				if not args.pieces:
					outdf = cluster_completeness(intype, cluster_map, inf2)
					outdf = outdf.round(decimals=3)
					outdf.to_csv(outf)
				else:
					print('\nOk, processing input file in pieces...\n')
					inkey = generate_index_list(inf2)
					c_list = list(cluster_map.keys())
					ct = len(c_list)
					print('Found %d clusters...' % ct)
					mat = np.zeros((ct, ct))  # initializes an array of the dimensions necessary to fit all cluster results
					j = 0
					for cluster in c_list:
						grab = pick_a_cluster(inkey, cluster)
						# print(grab)
						with open(args.input, 'r') as inf3:
							mat = big_cluster_completeness(grab, inkey, cluster, cluster_map, c_list, inf3, mat, j)
						# print(mat)
						j += 1
					print('File processing complete; writing output file...\n')
					outdf = pd.DataFrame(mat, dtype=float)
					outdf.columns = c_list  # names the columns (and index, next line) according to clusters in the order they were processed
					outdf.index = c_list
					outdf.sort_index(axis=0, inplace=True)
					outdf.sort_index(axis=1, inplace=True)
					outdf = outdf.round(decimals=3)
					outdf.to_csv(outf)

if __name__ == '__main__':
	main()
