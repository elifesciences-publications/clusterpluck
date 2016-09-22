#!/usr/bin/env Python

import argparse
import sys
import numpy as np
import pandas as pd
from collections import defaultdict


# The arg parser
def make_arg_parser():
	parser = argparse.ArgumentParser(description='Generates a binary OFU matrix (genomes vs BGCs) from clustered clusters')
	parser.add_argument('-i', '--input', help='Input file: a CSV of the BGCs and hclust clusters from R.', required=True)
	parser.add_argument('-c', '--completes', help='Fractional completeness matrix CSV file')
	parser.add_argument('-o', '--output', help='Where to save the output csv; default to screen', required=False, default='-')
	return parser


def outer(size):
	return lambda: np.zeros(size)

	dd['x'][7] += 1
	dd['x']


def main():
	parser = make_arg_parser()
	args = parser.parse_args()
	# Parse command line
	with open(args.input, 'r') as inf:
		hclus = pd.read_csv(inf, sep=',', header=0, index_col=0)
		size = hclus[1].max(0)
		outer(size)
		dd = defaultdict(outer(size))
		for line in hclus:
			strain = '_'.join(line.split('_')[0:2])
			if line.startswith(strain):
				pass
			else:
				pass



	df = pd.DataFrame(dd)

if __name__ == '__main__':
	main()
