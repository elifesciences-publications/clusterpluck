#!/usr/bin/env Python

import argparse
import os
import sys
import pandas as pd
from clusterpluck.scripts.cluster_dictionary import build_cluster_map
from clusterpluck.scripts.orfs_in_common import generate_index_list
from clusterpluck.scripts.orfs_in_common import pick_a_cluster


# The arg parser
def make_arg_parser():
	parser = argparse.ArgumentParser(description='Split a really big (or small) ORF matrix into pieces for reasonable processing. [-s] flag then compiles the results in a single directory')
	parser.add_argument('-i', '--input', help='Input is the ORF matrix CSV file.', default='-')
	parser.add_argument('-m', '--mpfa', help='The multi-protein fasta file (.mpfa) from which to build the dictionary')
	parser.add_argument('-b', '--bread', help='Where to find the cluster information in the header for the sequence (default="ref|,|")', default='ref|,|')
	parser.add_argument('-o', '--output', help='Where to save the output csv', required=False, default='cluster_chunk')
	parser.add_argument('-c', '--cuts', help='How many pieces the matrix should be split into', required=False, default=10)
	parser.add_argument('-s', '--synthesize', help='Run with this flag to put Humpty Dumpty back together (all csv in cwd!)', required=False, action='store_true', default=False)
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


def main():
	parser = make_arg_parser()
	args = parser.parse_args()
	if args.synthesize:
		final_df = synthesize_chunks()
		with open(args.output, 'w') as outf:
			final_df.to_csv(outf)
			print('\nMerged data written to file... exiting...\n')
			sys.exit()
	# Parse command line
	with open(args.mpfa, 'r') as inf:
		# Generates dictionary with each unique 'refseq_cluster' as keys, ORFs as values
		cluster_map = build_cluster_map(inf, bread=args.bread)
	with open(args.input, 'r') as in_csv:
		print('\nOk, processing input file in pieces...\n')
		inkey = generate_index_list(in_csv)
	c_list = list(cluster_map.keys())
	ct = len(c_list)
	n = int(args.cuts)
	print('Found %d clusters... Making %d cuts...' % (ct, n))
	# g = orfct / cuts
	# if not g.is_integer():
	# 	g = int(g) + 1
	bcl = [c_list[i:i + n] for i in range(0, len(c_list), n)]
	print('\nMaster list generated... now doing the splits!')
	p = 1
	for c in bcl:
		grab_chunk = []
		for cluster in list(c):
			grab = pick_a_cluster(inkey, cluster)  # uses the name of the cluster to get a list of all orfs for a particular unique cluster
			grab_chunk.extend(grab)
		with open(args.input, 'r') as inf3:
			mx = pd.read_csv(inf3, sep=',', header=0, usecols=grab_chunk, engine='c')  # loads in only the columns from the grab list, i.e. all cols for a unique cluster
		mx.index = inkey  # reindexes the df with the orf labels after importing specific columns with usecols
		# data_to_pool.append(mx)  # create the list of dfs to map over for multiprocessing
		outf = args.output
		if outf.endswith('.csv'):
			outf.replace('.csv', '')
		outf = '_'.join([args.output, str(p), '.csv'])
		mx.to_csv(outf)
		print('\nSaved a matrix chunk...')
		p += 1
		# if __name__ == '__main__':
			# print('\nSending data to Workers... work, Workers, work!')
			# results = list(futures.map(partial(parallel_clustermean, c_list=c_list), data_to_pool))
			# print('\nFile processing complete; writing output file...\n')
			# del data_to_pool
		# outdf.sort_index(axis=0, inplace=True)  # ensure that the clusters are in order on cols and rows
		# outdf.sort_index(axis=1, inplace=True)


if __name__ == '__main__':
	main()
