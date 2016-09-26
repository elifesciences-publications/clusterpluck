#!/usr/bin/env Python

import argparse
import sys
import pandas as pd
from scipy.cluster.hierarchy import linkage
from scipy.cluster.hierarchy import cut_tree


# The arg parser
def make_arg_parser():
	parser = argparse.ArgumentParser(description='Generates a binary OFU matrix (genomes vs BGCs) from clustered clusters')
	parser.add_argument('-i', '--input', help='Input file: A scores matrix of strain_cluster comparisons', required=True)
	parser.add_argument('-o', '--output', help='Where to save the output csv; default to screen', required=False, default='-')
	parser.add_argument('-h', '--height', help='At what height to cut the tree', required=False, default=500)
	return parser


# initialize the dict with zeros for each OFU
def process_hierarchy(inf, h):
	df = pd.read_csv(inf, header=0, index_col=0)
	df = df.fillna(0)
	strains = df.index
	li = linkage(df, method='ward', metric='euclidean')
	hclus = cut_tree(li, height=h)
	hclus = pd.DataFrame(li)
	hclus = hclus.set_index(strains)
	return hclus


def main():
	parser = make_arg_parser()
	args = parser.parse_args()
	# Parse command line

	with open(args.input, 'r') as inf:
		h = args.height
		hclus = process_hierarchy(inf, h)
		with open(args.output, 'w') if args.output != '-' else sys.stdout as outf:
			hclus.to_csv(outf)

if __name__ == '__main__':
	main()
