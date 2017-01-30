#!/usr/bin/env Python

import argparse
import sys
import os
import pandas as pd
from clusterpluck.wrappers.run_shogun import run_shogun
from clusterpluck.scripts.otu_x_ofu import match_tables
from clusterpluck.scripts.otu_x_ofu import multiply_tables
from multiprocessing import cpu_count


# The arg parser
def make_arg_parser():
	parser = argparse.ArgumentParser(description='Run SHOGUN lca with UTree, then compute 1-100% OFU tables.')
	parser.add_argument('-i', '--input', help='Directory containing the FASTA files (.fna) to run.', required=True)
	parser.add_argument('-p', '--threads', help='Number of threads to use in SHOGUN.', required=False, default=36)
	parser.add_argument('-u', '--utree', help='Path to the UTree database (.ctr file).', required=True)
	parser.add_argument('-c', '--cp_resources', help='Path to the ofu_keys directory, OFU files names as *_[%ID]ofu.csv', required=True)
	parser.add_argument('-m', '--multiples', help='When multiple strains for species-level match:'
												  '[average] the OFU tallies,'
												  '[summarize] all the possible OFUs,'
												  'select only [universal] OFUs,'
												  'select only OFUs in a [majority] of strains', required=False, default='majority')
	parser.add_argument('-o', '--output', help='Where to save the output data; default = cwd', required=False, default='shofun_utree_results')
	return parser


def main():
	parser = make_arg_parser()
	args = parser.parse_args()
	# Parse command line
	cpus = int(args.threads)
	cpus_avail = cpu_count()
	opt = args.multiples
	if cpus_avail > cpus:
		print('The number of processors requested is not available; maximum on this system is %s' % cpus_avail)
		sys.exit()
	indir = args.input
	if args.output:
		outdir = args.output
	else:
		outdir = os.path.join('shOFUn_utree')
	utree = args.utree
	cp_dir = args.cp_resources
	run_shogun(indir, outdir, utree, cpus)
	print('\nSHOGUN finished!... calculating OFU profiles...')
	for inofu in os.listdir(cp_dir):
		inofu_num = inofu.split('_')[-1]
		inofu_num = inofu_num.split('.')[0]
		with open(os.path.join(outdir, 'taxon_counts.csv'), 'r') as taxons:
			ofu_matched = match_tables(taxons, inofu, opt)
		with open(os.path.join(outdir, 'taxon_counts.csv'), 'r') as taxons:
			ofu_table = multiply_tables(taxons, ofu_matched)
			ofu_table = ofu_table.loc[:, (ofu_table != 0).any(axis=0)]  # removes ofus with all zeros
			print('Final ofu profile dimensions = ', ofu_table.shape[0], ',', ofu_table.shape[1], '\n')
		os.mkdir(os.path.join(outdir, 'ofu_profiles'))
		ofu_output = os.path.join(outdir, 'ofu_profiles')
		with open(os.path.join(ofu_output, 'profile_at_', inofu_num, '.csv'), 'w') as outf:
			ofu_table = ofu_table.round(decimals=2)
			ofu_table.to_csv(outf)
	print('\nOFU profiles written to file\n')


if __name__ == '__main__':
	main()