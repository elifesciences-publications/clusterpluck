#!/usr/bin/env Python

import argparse
import sys
from collections import defaultdict
from ninja_utils.parsers import FASTA
from ninja_utils.utils import find_between

# >ncbi_tid|206672|ref|NC_004307.2_cluster004_ctg1_orf02030|organism|Bifidobacterium_longum|
# The arg parser


def make_arg_parser():
	parser = argparse.ArgumentParser(
		description='Build a dictionary to store the list of ORFs for clusters in each genome')
	parser.add_argument('-i', '--input', help='Input is a multi-cluster protein FASTA file.', default='-')
	parser.add_argument('-b', '--bread', help='Where to find the header for the sequence (default="ref|,|")', default='ref|,|')
	return parser


# define the dictionary function
def build_cluster_map(inf, bread='ref|,|'):
	begin,end = bread.split(',')
	cluster_map = defaultdict(set)
	fasta_gen = FASTA(inf)
	for header, sequence in fasta_gen.read():
		if '.cluster' in header:
			header = header.replace('.cluster', '_cluster')
		ref = find_between(header, begin, end)
		header_split = ref.split('_')
		key = '_'.join(header_split[:3])
		value = '_'.join(header_split[-2:])
		cluster_map[key].add(value)
	return cluster_map


def main():
	parser = make_arg_parser()
	args = parser.parse_args()

	# parse command line
	with open(args.input, 'r') if args.input != '-' else sys.stdin as inf:
		cluster_map = build_cluster_map(inf, bread=args.bread)
		# debug
		# how many clusters there are
		print(len(list(cluster_map.keys())))
		print(list(cluster_map.keys()))
		# how many ORFs per cluster for all clusters
		print([len(value) for value in cluster_map.values()])
		# list all ORFs in a particular cluster
		# print(cluster_map['NC_010816.1_cluster001'])
		# how many ORFs in a particular cluster
		# print(len(cluster_map['NC_000000.0_cluster000]))
if __name__ == '__main__':
	main()
