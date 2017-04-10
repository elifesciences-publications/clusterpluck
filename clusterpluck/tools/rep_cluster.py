#!/usr/bin/env Python

import argparse
import sys
import os
import pandas as pd
from multiprocessing import cpu_count


# The arg parser
def make_arg_parser():
	parser = argparse.ArgumentParser(
		description='Builds a phylogenetic tree (using FastTree) from a representative OFU sequence '
		'(default = longest) at a specified ID level')
	parser.add_argument('-i', '--input',
						help='Input file: The blastp output file (b6 format)', required=True)
	parser.add_argument('--ofu',
						help='OFU txt file where the first column are the clusters matching those in the b6 input file', required=True)
	parser.add_argument('-c', '--cpus',
						help='Number of cpus to use', required=False)
	parser.add_argument('-u', '--underscore',
						help='For clustersuck, the underscore position (integer) on which to define a cluster', required=False, default=3)
	parser.add_argument('-q', '--quiet',
						help='Do not print all clustersuck output to screen', action='store_true', required=False, default=True)
	parser.add_argument('-o', '--output',
						help='Output txt file, otherwise write to screen', required=False, default='-')
	return parser


def rep_cluster_pick(in_b6, und, cpus, quiet):
	temp_result_m = os.path.join('temp_matrix.csv')
	if quiet:
		os.system(' '.join(['clustersuck', in_b6, temp_result_m, str(und), 'tempfilter.txt', str(cpus), '> /dev/null']))
	else:
		os.system(' '.join(['clustersuck', in_b6, temp_result_m, str(und), 'tempfilter.txt', str(cpus)]))
	rep_pick = pd.read_csv(temp_result_m, header=0, index_col=0)
	if rep_pick.shape[0] == 1:
		rep_pick = str(list(rep_pick.columns)[0])
	else:
		rep_pick = rep_pick.sum(axis=0)
		rep_pick = rep_pick.to_dict()
		rep_pick = max(rep_pick, key=rep_pick.get)
	os.remove(temp_result_m)
	return rep_pick


def main():
	parser = make_arg_parser()
	args = parser.parse_args()
	# Parse command line
	quiet = args.quiet
	und = args.underscore
	if args.cpus:
		cpus = int(args.cpus)
	else:
		cpus = cpu_count() / 2
	in_b6 = args.input
	with open(args.ofu, 'r') as ofu_infile:
		cluster_filter = pd.read_table(ofu_infile, header=0, index_col=0, sep='\t')
		cluster_filter = list(cluster_filter.index)
	with open('tempfilter.txt', 'w') as tempfilt:
		for n in cluster_filter:
			tempfilt.write(n + '\n')
	rep_pick = str(rep_cluster_pick(in_b6, und, cpus, quiet))
	os.remove('tempfilter.txt')
	with open(args.output, 'w') if args.output != '-' else sys.stdout as outf:
		outf.write(rep_pick)


if __name__ == '__main__':
	main()
