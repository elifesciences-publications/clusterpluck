#!/usr/bin/env Python

import argparse
import sys
import csv
import numpy as np
import pandas as pd
from collections import defaultdict


# The arg parser
def make_arg_parser():
	parser = argparse.ArgumentParser(description='Generates a binary OFU matrix (genomes vs BGCs) from clustered clusters')
	parser.add_argument('-i', '--input', help='Input file: a CSV of the BGCs and hclust clusters from R.', required=True)
	parser.add_argument('-o', '--output', help='Where to save the output csv; default to screen', required=False, default='-')
	return parser


# initialize the dict with zeros for each OFU
def outer(size):
	return lambda: np.zeros(size, dtype=np.int)


def cluster_ofus(inf2, dd):
	ofu_names = []
	hcsv = csv.reader(inf2, delimiter=',')
	next(hcsv, None)
	for line in hcsv:
		name = line[0]
		refseq_id = '_'.join(name.split('_')[0:2])
		ofu = int(line[1])
		ofu_name = ('ofu', str('%03d' % int(ofu)))
		ofu_names.append('_'.join(ofu_name))
		dd[refseq_id][ofu] += 1
	df = pd.DataFrame.from_dict(dd)
	df = df.T
	df.drop([0], axis=1, inplace=True)
	ofu_cols = list(set(ofu_names))
	ofu_cols.sort()
	df.columns = ofu_cols
	df.sort_index(axis=0)
	df.sort_index(axis=1)
	return df


def main():
	parser = make_arg_parser()
	args = parser.parse_args()
	# Parse command line
	with open(args.input, 'r') as inf:
		hclus = pd.read_csv(inf, sep=',', header=0, index_col=0)
		# print(hclus)
		size = hclus.max(0)[0]
		size += 1
		fill = outer(size)
		dd = defaultdict(fill)
		with open(args.input, 'r') as inf2:
			df = cluster_ofus(inf2, dd)
	with open(args.output, 'w') if args.output != '-' else sys.stdout as outf:
		df.to_csv(outf)


if __name__ == '__main__':
	main()
