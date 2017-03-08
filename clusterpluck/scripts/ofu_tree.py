#!/usr/bin/env Python

import argparse
import sys
import os
import shutil
import pandas as pd
from collections import defaultdict
from multiprocessing import cpu_count
from clusterpluck.tools.h_clustering import process_hierarchy
from clusterpluck.wrappers.run_makeblastdb import run_makeblastdb
from clusterpluck.wrappers.run_blastp import run_blastp
from ninja_utils.parsers import FASTA
from clusterpluck.tools.suppress_print import suppress_stdout

# The arg parser
def make_arg_parser():
	parser = argparse.ArgumentParser(
		description='Builds a phylogenetic tree (using FastTree) from a representative OFU sequence '
		'(default = longest) at a specified ID level')
	parser.add_argument('-i', '--input',
						help='Input file: A percent identity (pident) scores matrix of strain_cluster comparisons', required=True)
	parser.add_argument('-t', '--height',
						help='The similarity/identity level at which you want BGCs summarized within the tree (0-100)',
						required=True, default=70)
	parser.add_argument('-mpfa',
						help='The protein fasta resource containing cluster amino acid sequences for this ClusterPluck database', required=True)
	parser.add_argument('-o', '--output',
						help='Directory in which to save the cluster information files (default = cwd)', required=False, default='.')
	parser.add_argument('-no_cleanup',
						help='Keep all intermediate temp and log files', action='store_true', required=False, default=False)
	parser.add_argument('-p', '--cpus',
						help='Number of cpus to use - relevant for blastp and clusterparse steps', required=False)
	parser.add_argument('-u', '--underscore',
						help='For clusterparse, the underscore position (int) on which to define a cluster', required=False, default=4)
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

#
# def ofu_dnasequences(inf_d, bgc, ofu_name, dna_outf):
# 	fasta_gen = FASTA(inf_d)
# 	bgc = str(bgc)
# 	for header, sequence in fasta_gen.read():
# 		if '.cluster' in header:
# 			header = header.replace('.cluster', '_cluster')
# 		if bgc in header:
# 			dna_outf.write(''.join(['>', ofu_name, ' ', header, '|SPLIT_HERE|']))
# 			dna_outf.write(''.join([sequence, '\n']))
# 	return dna_outf


def compile_ofu_sequences(inf_m, bgc, ofu_name, aa_outf):
	mpfa_gen = FASTA(inf_m)
	bgc = str(bgc)
	for header, sequence in mpfa_gen.read():
		if '.cluster' in header:
			header = header.replace('.cluster', '_cluster')
		if bgc in header:
			aa_outf.write(''.join(['>', ofu_name, '_', header, '\n']))
			aa_outf.write(''.join([sequence, '\n']))
	return aa_outf


def rep_seqs(temppath, cut_h, outpath, tempdir, und, cpus):
	#with suppress_stdout():
	rep_setfile = os.path.join(temppath, ''.join(['repseqs_id', cut_h, '.mpfa']))
	rep_set = open(rep_setfile, 'w')
	for f in os.listdir(temppath):
		if f.endswith('aasequences.mpfa'):
			fp = os.path.join(temppath, f)
			# print(f)
			tempdb = ''.join([f.split('aaseq')[0], 'db'])
			tempdbp = os.path.join(temppath, tempdb)
			tempdbn = ''.join([tempdbp, '/tempdb'])
			temp_result = os.path.join(temppath, ''.join(['blast', f.split('aaseq')[0], '.b6']))
			temp_result_m = os.path.join(temppath, ''.join(['blast_mat_', f.split('aaseq')[0], '.csv']))
			# os.system(' '.join(['makeblastdb -in', fp, '-hash_index', '-out', tempdbn, '-dbtype prot']))
			run_makeblastdb(fp, tempdbn)
			#os.system(' '.join(['blastp -query', fp, '-out', temp_result, '-db', tempdbn, '-max_hsps 1 -outfmt 6 -evalue 0.05 -num_threads', str(cpus)]))
			run_blastp(fp, temp_result, tempdbn, cpus, evalue=0.05)
			os.system(' '.join(['clusterparse', temp_result, temp_result_m, str(und), str(cpus)]))
			rep_pick = pd.read_csv(temp_result_m, header=0, index_col=0)
			if rep_pick.shape[0] == 1:
				rep_pick = str(list(rep_pick.columns)[0])
			else:
				rep_pick = rep_pick.sum(axis=1)
				rep_pick = rep_pick.to_dict()
				rep_pick = max(rep_pick, key=rep_pick.get)
			with open(fp, 'r') as ofu_mpfa:
				mpfa_gen = FASTA(ofu_mpfa)
				for header, sequence in mpfa_gen.read():
					if rep_pick in header:
						rep_set.write('>' + header + '\n')
						rep_set.write(sequence + '\n')
	rep_set.close()
	return rep_setfile


def main():
	parser = make_arg_parser()
	args = parser.parse_args()
	# Parse command line
	if args.cpus:
		cpus = int(args.cpus)
	else:
		cpus = cpu_count()
	und = args.underscore
	outpath = args.output
	cut_h = str(args.height)
	h = float(args.height)
	tempdir = ''.join(['temp_id', cut_h])
	h = 1 - (h / 100)
	if not os.path.isdir(outpath):
		os.mkdir(outpath)
		os.mkdir(os.path.join(outpath, tempdir))
		if not os.path.isdir(outpath):
			print('\nError creating output directory; check given path and try again\n')
			sys.exit()
	with open(args.input, 'r') as inf:
		hclus = process_hierarchy(inf, h)
	bgc_dd, num_ofus = ofu_dictionary(hclus)
	for i in range(1, (num_ofus + 1)):
		ofu_name = ''.join(['ofu', ('%05d' % i)])
		ofu_digits = '%05d' % i
		ofu_aaseqfile = ''.join([ofu_name, '_id', cut_h, '_aasequences.mpfa'])
		aa_outf = open(os.path.join(outpath, tempdir, ofu_aaseqfile), 'w')
		bgcs = bgc_dd[ofu_digits]
		for bgc in bgcs:
			with open(args.mpfa, 'r') as inf_m:
				aa_outf = compile_ofu_sequences(inf_m, bgc, ofu_name, aa_outf)
		aa_outf.close()
	print('Finished compiling OFU DNA sequences... Now picking representative sequences...\n')
	temppath = os.path.join(outpath, tempdir)
	rep_setfile = rep_seqs(temppath, cut_h, outpath, tempdir, und, cpus)
	repdb = ''.join(['repseq_blastdb'])
	repdbp = os.path.join(outpath, tempdir, repdb)
	repdbn = ''.join([repdbp, 'repdb'])
	rep_result = os.path.join(outpath, tempdir, 'repseq_blast.b6')
	rep_result_m = os.path.join(outpath, tempdir, 'repset_matrix.csv')
	# with suppress_stdout():
	os.system(' '.join(['makeblastdb -in', rep_setfile, '-hash_index', '-out', repdbn, '-dbtype prot']))
	os.system(' '.join(['blastp -query', rep_setfile, '-out', rep_result, '-db', repdbn, '-max_hsps 1 -outfmt 6 -evalue 0.05 -num_threads', str(cpus)]))
	os.system(' '.join(['clusterparse', rep_result, rep_result_m, str(und), str(cpus)]))
	newick_tree = os.path.join(outpath, ''.join(['OFU_tree_id', cut_h, '.tree']))
	# print(' '.join(['make_ofu_tree.R', rep_result_m, newick_tree]))
	os.system(' '.join(['make_ofu_tree.R', rep_result_m, newick_tree]))
	if not args.no_cleanup:
		print('Alrighty, cleaning up temp files...\n')
		shutil.rmtree((os.path.join(outpath, tempdir)), ignore_errors=False, onerror=None)


if __name__ == '__main__':
	main()
