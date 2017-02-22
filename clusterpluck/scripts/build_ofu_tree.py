#!/usr/bin/env Python

import argparse
import sys
import os
import shutil
from collections import defaultdict
from clusterpluck.tools.h_clustering import process_hierarchy
from ninja_utils.parsers import FASTA

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
	parser.add_argument('-dna_fasta',
						help='The multi-fasta resource containing cluster DNA sequences for this ClusterPluck database', required=True)
	parser.add_argument('-o', '--output',
						help='Directory in which to save the cluster information files (default = cwd)', required=False, default='.')
	parser.add_argument('-no_cleanup',
						help='Keep all intermediate temp and log files', action='store_true', required=False, default=False)
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


def ofu_dnasequences(inf_d, bgc, ofu_name, dna_outf):
	fasta_gen = FASTA(inf_d)
	bgc = str(bgc)
	for header, sequence in fasta_gen.read():
		if '.cluster' in header:
			header = header.replace('.cluster', '_cluster')
		if bgc in header:
			dna_outf.write(''.join(['>', ofu_name, ' ', header, '|SPLIT_HERE|']))
			dna_outf.write(''.join([sequence, '\n']))
	return dna_outf


def rep_seqs(inpath, cut_h):
	align_out = os.path.join(inpath, ''.join(['repseqs_id', cut_h, '.fna']))
	to_align = open(align_out, 'w')
	for f in os.listdir(inpath):
		if f.endswith('dnasequences.fna'):
			f = os.path.join(inpath, f)
			print(f)
			longrep = max(open(f, 'r'), key=len)
			to_align.write(longrep.replace('|SPLIT_HERE|', '\n'))
	to_align.close()
	return align_out


def main():
	parser = make_arg_parser()
	args = parser.parse_args()
	# Parse command line
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
		ofu_dnaseqfile = ''.join([ofu_name, '_id', cut_h, '_dnasequences.fna'])
		dna_outf = open(os.path.join(outpath, tempdir, ofu_dnaseqfile), 'w')
		bgcs = bgc_dd[ofu_digits]
		for bgc in bgcs:
			with open(args.dna_fasta, 'r') as inf_d:
				dna_outf = ofu_dnasequences(inf_d, bgc, ofu_name, dna_outf)
		dna_outf.close()
	print('Finished compiling OFU DNA sequences... Now picking representative sequences...\n')
	inpath = os.path.join(outpath, tempdir)
	to_align_out = rep_seqs(inpath, cut_h)
	print('Running alignment with MUSCLE...\n')
	aligned = os.path.join(outpath, ''.join(['OFU_alignment_id', cut_h, '.fa']))
	alignment_log = os.path.join(outpath, tempdir, ''.join(['OFU_alignment_id', cut_h, 'log.log']))
	os.system(' '.join(['muscle', '-in', to_align_out, '-out', aligned, '-log', alignment_log, '-maxiters 3', '-diags1']))
	# DEBUG
	# print(' '.join(['muscle', '-in', to_align_out, '-out', aligned, '-log', alignment_log]))
	print('Alignment completed, building tree...\n')
	newick_tree = os.path.join(outpath, ''.join(['OFU_tree_id', cut_h, '.tree']))
	tree_log = os.path.join(outpath, tempdir, ''.join(['OFU_tree_id', cut_h, 'log.log']))
	os.system(' '.join(['FastTree', '-log', tree_log, '-nt', '-gtr', '-pseudo', aligned, '>', newick_tree]))
	# DEBUG
	# print(' '.join(['FastTree', '-log', tree_log, '-nt', '-gtr', '-pseudo', aligned, '>', newick_tree]))
	print('Tree written to Newick Format successfully!')
	if not args.no_cleanup:
		print('Alrighty, cleaning up temp files...\n')
		shutil.rmtree((os.path.join(outpath, tempdir)), ignore_errors=False, onerror=None)


if __name__ == '__main__':
	main()
