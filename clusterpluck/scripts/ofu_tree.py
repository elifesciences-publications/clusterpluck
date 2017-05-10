#!/usr/bin/env Python

import argparse
import sys
import os
import shutil
import pandas as pd
from collections import defaultdict
from itertools import repeat
from multiprocessing import Pool
from multiprocessing import cpu_count
from clusterpluck.tools.h_clustering import process_hierarchy
from clusterpluck.tools.product_profiler import type_dict


# The arg parser
def make_arg_parser():
	parser = argparse.ArgumentParser(
		description='Builds a phylogenetic tree (using FastTree) from a representative OFU sequence '
		'(default = longest) at a specified ID level')
	parser.add_argument('-i', '--input',
						help='Input file: The blastp output file (b6 format)', required=True)
	parser.add_argument('-s', '--scores',
						help='Optional input file: The pre-processed (with clustersuck) scores matrix', required=False)
	parser.add_argument('-t', '--height',
						help='The similarity/identity level at which you want BGCs summarized within the tree (0-100)',
						required=True, default=70)
	parser.add_argument('-m', '--method',
						help='What clustering method to use on the distance matrix: single, complete, average, weighted, centroid, median, ward.',
						required=False, default='average')
	parser.add_argument('-o', '--output',
						help='Directory in which to save the cluster information files (default = cwd)', required=False, default='.')
	parser.add_argument('-no_cleanup',
						help='Keep all intermediate temp and log files', action='store_true', required=False, default=False)
	parser.add_argument('-c', '--cpus',
						help='Number of cpus to use', required=False)
	parser.add_argument('-u', '--underscore',
						help='For clustersuck, the underscore position (integer) on which to define a cluster', required=False, default=4)
	parser.add_argument('-quiet',
						help='Do not print all clustersuck output to screen', action='store_true', required=False, default=False)
	parser.add_argument('-y', '--types',
						help='If an OFU-to-representative-cluster-product-type table is desired, provide the path to the strain product types csv.', required=False)
	return parser


# Create the OFU dictionary from the clustered results
def ofu_dictionary(hclus):
	bgc_dd = defaultdict(list)
	for value, key in hclus.itertuples(index=True):
		key = str('%05d' % key)
		bgc_dd[key].extend(list([value]))
	num_ofus = len(bgc_dd)
	print('\nPreparing to build tree for %d OFUs...\n' % num_ofus)
	return bgc_dd, num_ofus


def rep_cluster_pick(rep_ofu_file, args_list):
	in_b6 = args_list[0]
	temppath = args_list[1]
	quiet = args_list[2]
	und = args_list[3]
	cpus_clus = args_list[4]
	rep_file_path = os.path.join(temppath, rep_ofu_file)
	temp_result_m = os.path.join(temppath, ''.join([rep_ofu_file.split('_filter')[0], '_matrix.csv']))
	if quiet:
		os.system(' '.join(['clustersuck', in_b6, temp_result_m, str(und), rep_file_path, str(cpus_clus), '> /dev/null']))
	else:
		os.system(' '.join(['clustersuck', in_b6, temp_result_m, str(und), rep_file_path, str(cpus_clus)]))
	rep_pick = pd.read_csv(temp_result_m, header=0, index_col=0)
	if rep_pick.shape[0] == 1:
		rep_pick = str(list(rep_pick.columns)[0])
	else:
		rep_pick = rep_pick.sum(axis=0)
		rep_pick = rep_pick.to_dict()
		rep_pick = max(rep_pick, key=rep_pick.get)
	os.remove(temp_result_m)
	return rep_pick


def relabeler(rep_matrix, bgc_dd):
	rep_df = pd.read_csv(rep_matrix, header=0, index_col=0)
	rep_df.sort_index(axis=0, inplace=True)
	rep_df.sort_index(axis=1, inplace=True)
	names = list(rep_df.columns)
	# print(names[0:11])
	ofu_nums = []
	for n in names:
		ofu_nums.extend([key for key, value in bgc_dd.items() if n in value])
	ofu_ids = []
	for ofu in ofu_nums:
		ofu = int(ofu)
		full_id = str('%05d' % ofu)
		full_id = ''.join(['ofu', full_id])
		ofu_ids.append(full_id)
	# print(len(names))
	# print(len(ofu_nums))
	# print(len(ofu_ids))
	# print(ofu_ids[1:10])
	if not len(names) == len(ofu_ids):
		print('Error renaming data with OFUs; check for duplicate entries in data set')
		sys.exit()
	rep_df.columns = ofu_ids
	rep_df.index = ofu_ids
	return rep_df


def main():
	parser = make_arg_parser()
	args = parser.parse_args()
	# Parse command line
	if args.cpus:
		cpus = int(args.cpus)
	else:
		cpus = cpu_count()
	# Split the cpus rationally among the processes. At higher cpu numbers, ultimately I/O will reign.
	if cpus > 8:
		cpus_py = 4
		cpus_clus = int((cpus - 4) / 4)
	else:
		cpus_py = cpus
		cpus_clus = 1
	method = args.method
	in_b6 = args.input
	und = args.underscore
	outpath = args.output
	cut_h = str(args.height)
	h = float(args.height)
	h = 1 - (h / 100)
	# Make a temp location to store all the intermediate files
	tempdir = ''.join(['temp_id', cut_h])
	temppath = os.path.join(outpath, tempdir)
	if not os.path.isdir(outpath):
		os.mkdir(outpath)
	if not os.path.isdir(temppath):
		os.mkdir(os.path.join(outpath, tempdir))
	if not os.path.isdir(outpath):
		print('\nError creating output directory; check given path and try again\n')
		sys.exit()
	if args.scores:
		with open(args.scores, 'r') as inf:
			hclus = process_hierarchy(inf, h, method)
	# Make the scores matrix from blast result if none provided
	else:
		full_mat = os.path.join(temppath, 'full_scores_matrix.csv')
		os.system(' '.join(['clustersuck', in_b6, full_mat, str(und), 'NOFILTER', str(cpus)]))
		# Run hclus on the scores matrix, generating the OFUs themselves
		with open(full_mat, 'r') as inf:
			hclus = process_hierarchy(inf, h, method)
	bgc_dd, num_ofus = ofu_dictionary(hclus)
	for i in range(1, (num_ofus + 1)):
		ofu_name = ''.join(['ofu', ('%05d' % i)])
		ofu_digits = '%05d' % i
		ofu_clusterfile = ''.join([ofu_name, '_id', cut_h, '_filter.txt'])
		ofu_filter = open(os.path.join(temppath, ofu_clusterfile), 'w')
		bgcs = bgc_dd[ofu_digits]
		for bgc in bgcs:
			ofu_filter.write(bgc + '\n')
		ofu_filter.close()
	print('Finished compiling OFU clusters... Now picking best representatives...\n')
	rep_clusters = os.path.join(temppath, 'rep_clusters_%s_id.txt' % cut_h)
	ofu_rep_list = [f for f in os.listdir(temppath) if f.endswith('filter.txt')]
	quiet = args.quiet
	args_list = [in_b6, temppath, quiet, und, cpus_clus]
	# Run the representative cluster function in parallel
	with Pool(processes=cpus_py) as pool:
		results = pool.starmap(rep_cluster_pick, zip(ofu_rep_list, repeat(args_list)))
		pool.close()
		pool.join()
	with open(rep_clusters, 'w') as outf:
		for item in results:
			outf.write('%s\n' % item)
	rep_result_m = os.path.join(temppath, 'repset_matrix_%s.csv' % cut_h)
	# Run clustersuck on the representative clusters, to generate the representative scores matrix at chosen height
	os.system(' '.join(['clustersuck', in_b6, rep_result_m, str(und), rep_clusters, str(cpus)]))
	if args.types:
		type_dd, types = type_dict(args.types)
		ofu_types_fp = os.path.join(outpath, ''.join(['ofu_', cut_h, '_cluster_types.csv']))
		ofu_types = open(ofu_types_fp, 'w')
		ofu_types.write(',cluster_type\n')
		with open(rep_result_m, 'r') as rep_matrix:
			rep_df = pd.read_csv(rep_matrix, header=0, index_col=0)
		rep_df.sort_index(axis=0, inplace=True)
		rep_df.sort_index(axis=1, inplace=True)
		names = list(rep_df.columns)
		for n in names:
			p_type = str(type_dd[n])
			p_ofu = [key for key, value in bgc_dd.items() if n in value]
			print(p_ofu)
			p_ofu = int(p_ofu[0])
			full_ofu = str('%05d' % p_ofu)
			full_ofu = ''.join(['ofu', full_ofu])
			ofu_types.write(full_ofu + ',' + p_type + '\n')
		ofu_types.close()
	with open(rep_result_m, 'r') as rep_matrix:
		rep_df = relabeler(rep_matrix, bgc_dd)
	rep_ofu_result = os.path.join(outpath, 'ofu_repset_matrix.csv')
	with open(rep_ofu_result, 'w') as ofu_outf:
		rep_df.to_csv(ofu_outf)
	newick_tree = os.path.join(outpath, ''.join(['OFU_tree_id', cut_h, '.tree']))
	# print(' '.join(['make_ofu_tree.R', rep_result_m, newick_tree]))
	# Make the tree with R from the representative cluster scores matrix
	os.system(' '.join(['make_ofu_tree.R', rep_ofu_result, newick_tree, method]))
	if not args.no_cleanup:
		print('\nAlrighty, cleaning up temp files...\n')
		shutil.rmtree((os.path.join(outpath, tempdir)), ignore_errors=False, onerror=None)


if __name__ == '__main__':
	main()
