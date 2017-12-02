#!/usr/bin/env Python

import argparse
import os
import sys
import numpy as np
from collections import defaultdict
from multiprocessing import cpu_count

import pandas as pd
from ninja_utils.parsers import FASTA
from clusterpluck.tools.suppress_print import suppress_stdout


# Arg parser
def make_arg_parser():
	parser = argparse.ArgumentParser(description='Given the processed shotgun reads, measure coverage of predicted BGCs in metagenome')
	parser.add_argument('-s', '--seqs', help="The combined seqs file (fasta) containing QC'd shotgun reads", required=True)
	parser.add_argument('--bgc_fna', help="The BGC DNA fasta file (all BGCs in one multi-fasta file)", required=False)
	parser.add_argument('--bgc_edx', help="The pre-made BURST database file)", required=False)
	parser.add_argument('--bgc_acx', help="The pre-made BURST accelerator file)", required=False)
	# parser.add_argument('-b', '--bread', help='Where to find the header for the sequence (default="ref|,|")', default='ref|,|')
	parser.add_argument('-o', '--output', help='Directory in which to save the results (default = cwd)', required=False, default='.')
	parser.add_argument('-m', '--metrics', help='Coverage metrics to report (default = all)', required=False, default='-')
	parser.add_argument('-c', '--ofu', help='Comma-separated list of the ofus (e.g. ofu00001,ofu00003) on which to provide information', required=False, type=str)
	parser.add_argument('--nt_cat', help='Path to the nt catalog, required for genbank ID clusters (i.e. antismash DB)', required=False)
	parser.add_argument('-t', '--threshold', help='Coverage threshold at which to accept a cluster alignment (default = 75%)', required=False, default=75, type=float)
	parser.add_argument('-p', '--threads', help='Number of processors to use for alignment (default = all available)', required=False, default='-', type=int)
	return parser


# def compile_ofu_dnasequences(inf_d, bgc, dna_outf):
# 	fasta_gen = FASTA(inf_d)
# 	bgc = str(bgc)
# 	for header, sequence in fasta_gen.read():
# 		if '.cluster' in header:
# 			header = header.replace('.cluster', '_cluster')
# 		if bgc in header:
# 			dna_outf.write(''.join(['>', header, '\n']))
# 			dna_outf.write(''.join([sequence, '\n']))
# 	return dna_outf


# Make a database if necessary
def create_bgc_db(bgc_db, temp_path, burst):
	acc = os.path.join(temp_path, 'bgc_db.acc')
	edb = os.path.join(temp_path, 'bgc_db.edb')
	os.system(' '.join([burst, '-r', bgc_db, '-a', acc, '-o', edb, '-d DNA -s']))
	os.system(' '.join([burst, '-a', acc, '-o /dev/null']))
	os.system(' '.join([burst, '-r', edb, '-o /dev/null']))
	os.remove(acc)
	os.remove(edb)
	acx = os.path.join(temp_path, 'bgc_db.acx')
	edx = os.path.join(temp_path, 'bgc_db.edx')
	if not os.path.isfile(acx) or not os.path.isfile(edx):
		raise ValueError('BURST databases not successfully generated: check input fasta and BURST installation')
	return acx, edx


# Run the BURST alignment
def align_bgcs(outpath, seqs, acx, edx, cpus):
	alignment = os.path.join(outpath, 'bgc_alignments.b6')
	if acx != '-':
		os.system(' '.join(['burst15 -m ALLPATHS -fr -hr -q', seqs, '-n -r', edx, '-o', alignment, '-t', cpus, '-i 0.9', '-a', acx]))
	else:
		os.system(' '.join(['burst15 -m ALLPATHS -fr -hr -q', seqs, '-n -r', edx, '-o', alignment, '-t', cpus, '-i 0.9']))
	print("BURST alignment completed")
	return alignment


def coverage_analysis(alignment, outpath, output_metrics):
	alignment_t = pd.read_table(alignment, header=None, index_col=0)
	# TODO: Coverage... or just use shogun for the alignment and coverage?
	np.array()
	coverage_result = os.path.join(outpath, 'coverage_result.txt')
	# Run the allpaths coverage parsing script with desired output parameters
	os.system(' '.join(['allpaths_are_belong_to_us', alignment, coverage_result, output_metrics]))
	return coverage_result


def main():
	parser = make_arg_parser()
	args = parser.parse_args()
	outpath = args.output
	seqs = args.seqs
	if args.threads == '-':
		cpus = cpu_count()
	else:
		cpus = int(args.threads)
	if not args.metrics == '-':
		output_metrics = str(args.metrics)
	else:
		args.metrics = 'all'
	if not os.path.isfile(seqs):
		raise ValueError('Sequences file not found')
	threshold = float(args.threshold)
	if not os.path.isdir(outpath):
		os.mkdir(outpath)
		if not os.path.isdir(outpath):
			raise ValueError('Error creating output directory; check given path and try again')
	temp_path = os.path.join(outpath, 'tempz')
	os.mkdir(temp_path)
	if args.bgc_fna:
		fna_bytes = int(os.path.getsize(args.bgc_fna))
		if fna_bytes < 1000000000:
			burst = 'burst12'
		else:
			burst = 'burst15'
		bgc_db = args.bgc_fna
		acx, edx = create_bgc_db(bgc_db, temp_path, burst)
	elif args.bgc_edx:
		edx = args.bgc_edx
		if args.bgc_acx:
			acx = args.bgc_acx
		else:
			acx = '-'
	elif not args.bgc_fna and not args.bgc_edx:
		raise ValueError('No reference sequences or database provided: please include something to align to')

	# Generate the alignment between sequences and BGCs
	alignment = align_bgcs(outpath, seqs, acx, edx, cpus)

	# TODO: Filter the alignment by predicted profiles?
	# Or just see what is there?

	coverage_result = coverage_analysis(alignment, outpath, output_metrics)


if __name__ == '__main__':
	main()
