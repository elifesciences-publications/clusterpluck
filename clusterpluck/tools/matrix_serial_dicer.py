#!/usr/bin/env Python

import argparse
import os
import sys
import pandas as pd
from clusterpluck.scripts.cluster_dictionary import build_cluster_map
from clusterpluck.scripts.orfs_in_common import generate_index_list
from clusterpluck.scripts.orfs_in_common import pick_a_cluster
from itertools import repeat
from multiprocessing import Pool
from multiprocessing import cpu_count


# The arg parser
def make_arg_parser():
	parser = argparse.ArgumentParser(description='Split a really big (or small) ORF matrix into pieces for reasonable processing. [-s] flag then compiles the results in a single directory')
	parser.add_argument('-i', '--input', help='Input is the ORF matrix CSV file.', default='-')
	parser.add_argument('-m', '--mpfa', help='The multi-protein fasta file (.mpfa) from which to build the dictionary')
	parser.add_argument('-b', '--bread', help='Where to find the cluster information in the header for the sequence (default="ref|,|")', default='ref|,|')
	parser.add_argument('-o', '--output', help='Where to save the output csv', required=False, default='cluster_chunk')
	parser.add_argument('-c', '--cutsize', help='Define the number of clusters per chunk', required=False, default=10)
	parser.add_argument('-s', '--synthesize', help='Run with this flag to put Humpty Dumpty back together (all csv in cwd!)', required=False, action='store_true', default=False)
	parser.add_argument('-p', '--cpus', help='Number of processors to use', required=False, default=4, type=int)
	return parser


def synthesize_chunks():
	filelist = os.listdir('.')
	filechunks = [n for n in filelist if n.endswith('.csv')]
	data_list = []
	for f in filechunks:
		with open(f, 'r') as inf:
			c = pd.read_csv(inf, header=0, index_col=0, engine='c')
			c.sort_index(axis=0, inplace=True)
			data_list.append(c)
	final_df = pd.concat(data_list, axis=1)
	final_df.sort_index(axis=1, inplace=True)
	return final_df
#
#
# def parallel_chunker(c, arglist):
# 	infile = arglist[0]
# 	inkey = arglist[1]
# 	grab_chunk = []
# 	for cluster in list(c):
# 		grab = pick_a_cluster(inkey, cluster)  # uses the name of the cluster to get a list of all orfs for a particular unique cluster
# 		grab_chunk.extend(grab)
# 	with open(infile, 'r') as inf3:
# 		mx = pd.read_csv(inf3, sep=',', header=0, usecols=grab_chunk, engine='c')  # loads in only the columns from the grab list, i.e. all cols for a unique cluster
# 		mx.index = inkey  # reindexes the df with the orf labels after importing specific columns with usecols
# 	return mx


def main():
	parser = make_arg_parser()
	args = parser.parse_args()
	# Parse command line
	if args.synthesize:
		final_df = synthesize_chunks()
		with open(args.output, 'w') as outf:
			final_df.to_csv(outf)
			print('\nMerged data written to file... exiting...\n')
			sys.exit()
	with open(args.mpfa, 'r') as inf:
		# Generates dictionary with each unique 'refseq_cluster' as keys, ORFs as values
		cluster_map = build_cluster_map(inf, bread=args.bread)
	with open(args.input, 'r') as in_csv:
		print('\nOk, processing input file...\n')
		big_df = pd.read_csv(in_csv, sep=',', header=0, index_col=0, engine='c')
	inkey = list(big_df.index)
		# inkey = generate_index_list(in_csv)
	c_list = list(cluster_map.keys())
	ct = len(c_list)
	n = int(args.cutsize)
	print('Found %d clusters... Making groups of %d clusters...' % (ct, n))
	# Make a list of lists of clusters, to guide the breaking up of the csv
	bcl = [c_list[i:i + n] for i in range(0, len(c_list), n)]
	print('\nMaster list generated... now doing the splits!')
	p = 1
	for c in bcl:
		grab_chunk = []
		for cluster in list(c):
			grab = pick_a_cluster(inkey, cluster)  # uses the name of the cluster to get a list of all orfs for a particular unique cluster
			grab_chunk.extend(grab)
		chunk_df = big_df[grab_chunk]
		# with open(args.input, 'r') as inf3:
		# 	mx = pd.read_csv(inf3, sep=',', header=0, usecols=grab_chunk, engine='c')  # loads in only the columns from the grab list, i.e. all cols for a unique cluster
		# mx.index = inkey  # reindexes the df with the orf labels after importing specific columns with usecols
		# data_to_pool.append(mx)  # create the list of dfs to map over for multiprocessing
		outf = args.output
		if outf.endswith('.csv'):
			outf.replace('.csv', '')
		outf = '_'.join([outf, str(p), '.csv'])
		chunk_df.to_csv(outf)
		print('\nSaved matrix chunk %d...' % p)
		p += 1
	# outdf.sort_index(axis=0, inplace=True)  # ensure that the clusters are in order on cols and rows
	# outdf.sort_index(axis=1, inplace=True)

if __name__ == '__main__':
	main()
