#!/usr/bin/env python
#
# Robin Shields-Cutler
# August 2016
# takes standard blast output TSV (outfmt 6-- *.b6 or *.txt, etc) and stores entries in dictionary,
# then writes to dataframe and exports as CSV

usage = 'blastp_to_matrix.py -i BLASTOUT.b6 -s SCORE_METHOD -t THRESHOLD -o OUTFILE.csv'

import argparse
import os
import csv
import pandas as pd
import numpy as np
import re
from collections import defaultdict


def make_arg_parser():
	parser = argparse.ArgumentParser(description='Convert blastp output txt table to a scores matrix in csv format')
	parser.add_argument('-i', '--input', help='The blast output file to process.', required=True, type=str)
	parser.add_argument('-s', '--score', help='Which score to enter into matrix: "pident", "evalue", or "bitscore"', required=False, type=str, default='bitscore')
	parser.add_argument('-t', '--threshold', help='The threshold (float) for entry into matrix.', required=False, type=float, default=1)
	parser.add_argument('-o', '--output', help='Where to put the output (CSV or h5)', required=False, type=str, default='blastp_matrixform.csv')
	parser.add_argument('-n', '--normalize', help='Normalize bitscore to score of self-self for each cluster (as 100).', action='store_true', required=False, default=False)
	return parser


def main():
	parser = make_arg_parser()
	args = parser.parse_args()

	sparse_blast_id_dict = defaultdict(dict)
	with open(args.input) as blast_inf:
		# next(blast_inf)
		blast_tsv = csv.reader(blast_inf, delimiter='\t')
	# line[0] query name, line[1] = reference name, line[2] = % match, line[10] = e-value, line[11] = bitscore
		if args.score == 'bitscore':
			self_match_dict = defaultdict(dict)
			for line in blast_tsv:
				p = re.compile(r'(\w+_[\w+\d+]*\.\d)(_\w+\d\d\d)(_ctg\d_orf\d+)')
				m = p.search(line[0])
				n = p.search(line[1])
				# mref = m.group(1, 2)
				# nref = n.group(1, 2)
				cname = ''.join(m.group(1, 2, 3))
				rname = ''.join(n.group(1, 2, 3))
				# # print(mref)
				bvalue = np.float(line[11])
				if cname == rname:
					self_match_dict[line[0]] = bvalue
				if bvalue > args.threshold:
					sparse_blast_id_dict[line[0]][line[1]] = bvalue
		# TODO: use the evalue of perfect matches to normalize the data
		elif args.score == 'evalue':
			for line in blast_tsv:
				# p = re.compile(r'(\w+_[\w+\d+]*\.\d)(_\w+\d\d\d)(_ctg\d_orf\d+)')
				# m = p.search(line[0])
				# n = p.search(line[1])
				# mref = m.group(1)
				# nref = n.group(1)
				# cname = ''.join(m.group(1, 2, 3))
				# rname = ''.join(n.group(1, 2, 3))
				# if mref == nref:
				# 	pass
				# else:
				evalue = np.float(line[10])
				if evalue < args.threshold:
					sparse_blast_id_dict[line[0]][line[1]] = evalue
		elif args.score == 'pident':
			for line in blast_tsv:
				# p = re.compile(r'(\w+_[\w+\d+]*\.\d)(_\w+\d\d\d)(_ctg\d_orf\d+)')
				# m = p.search(line[0])
				# n = p.search(line[1])
				# mref = m.group(1, 2)
				# nref = n.group(1, 2)
				# cname = ''.join(m.group(1, 2, 3))
				# rname = ''.join(n.group(1, 2, 3))
				# if mref == nref:
				# 	pass
				# else:
				ivalue = np.float(line[2])
				if ivalue > args.threshold:
					sparse_blast_id_dict[line[0]][line[1]] = ivalue

	df = pd.DataFrame.from_dict(sparse_blast_id_dict)
	df.sort_index(axis=0, inplace=True)
	df.sort_index(axis=1, inplace=True)
	# print(df.shape[0])
	if args.normalize:
		vals = []
		for cluster in list(df.columns):
			vals.append(self_match_dict[cluster])
		df = df / vals * 100
		# print(len(vals))
	# Check if a matrix is symmetric
	# arr = df.values
	# print((arr.transpose() == -arr).all())
	if args.output.endswith('.csv'):
		df.to_csv(args.output)
	else:
		df.to_hdf(args.output, 'table')


if __name__ == '__main__':
	main()
