#!/usr/bin/env Python

import argparse
import sys
import pandas as pd


# The arg parser
def make_arg_parser():
	parser = argparse.ArgumentParser(description='multiplies each cluster x cluster average score by the fractional completeness')
	parser.add_argument('-s', '--scores', help='Average scores matrix CSV file.', required=True)
	parser.add_argument('-c', '--completes', help='Fractional completeness matrix CSV file')
	parser.add_argument('-o', '--output', help='Where to save the output csv; default to screen', required=False, default='-')
	return parser


def clustered(scores_inf, completes_inf):
	scores_df = pd.read_csv(scores_inf, sep=',', header=0, index_col=0)
	completes_df = pd.read_csv(completes_inf, sep=',', header=0, index_col=0)
	if scores_df.shape[0] != completes_df.shape[0] or scores_df.shape[1] != completes_df.shape[1]:
		print('\nWarning, input matrices do not have the same dimensions; unmatched data filled with NaN.\n')
	else:
		print('\nMatrices have equal dimensions... performing multiplication...\n')
	clustered_mx = scores_df.multiply(completes_df)
	clustered_mx = clustered_mx.round(decimals=1)
	clustered_mx.sort_index(axis=0)
	clustered_mx.sort_index(axis=1)
	# print(clustered_mx)
	return clustered_mx


def main():
	parser = make_arg_parser()
	args = parser.parse_args()
	# Parse command line
	with open(args.scores, 'r') as scores_inf:
		with open(args.completes, 'r') as completes_inf:
			with open(args.output, 'w') if args.output != '-' else sys.stdout as outf:
				clustered_mx = clustered(scores_inf, completes_inf)
				# Write the output csv file, if defined
				clustered_mx.to_csv(outf)

if __name__ == '__main__':
	main()
