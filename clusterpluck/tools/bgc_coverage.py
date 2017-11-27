#!/usr/bin/env Python

import argparse
import os
import sys
from collections import defaultdict

import pandas as pd
from ninja_utils.parsers import FASTA
from clusterpluck.tools.suppress_print import suppress_stdout


# The arg parser
def make_arg_parser():
	parser = argparse.ArgumentParser(description='Given the processed shotgun reads, measure coverage of predicted BGCs in metagenome')
	parser.add_argument('-s', '--seqs', help="The combined seqs file (fasta) containing QC'd shotgun reads", required=True)
	parser.add_argument('-d', '--bgc_fna', help="The BGC DNA fasta file (each BGC linearized then all BGCs in one multi-fasta file)", required=True)
	parser.add_argument('-b', '--bread', help='Where to find the header for the sequence (default="ref|,|")', default='ref|,|')
	parser.add_argument('-o', '--output', help='Directory in which to save the results (default = cwd)', required=False, default='.')
	parser.add_argument('-c', '--ofu', help='Comma-separated list of the ofus (e.g. ofu00001,ofu00003) on which to provide information', required=False, type=str)
	parser.add_argument('--nt_cat', help='Path to the nt catalog, required for genbank ID clusters (i.e. antismash DB)', required=False, default='-')
	parser.add_argument('-t', '--threshold', help='Coverage threshold at which to accept a cluster alignment (default = 75%)', required=False, default=75, type=float)
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


def create_bgc_db(bgc_db, temp_path):
	acc = os.path.join(temp_path, 'bgc_db.acc')
	edb = os.path.join(temp_path, 'bgc_db.edb')
	os.system(' '.join(['burst -r', bgc_db, '-a', acc, '-o', edb, '-f -d -s']))
	os.system(' '.join(['burst -r', edb, '-a', acc, '-o /dev/null']))
	os.system(' '.join(['burst -r', edb, '-o /dev/null']))
	os.remove(acc)
	os.remove(edb)
	acx = os.path.join(temp_path, 'bgc_db.acx')
	edx = os.path.join(temp_path, 'bgc_db.edx')
	if not os.path.isfile(acx) or not os.path.isfile(edx):
		raise ValueError('BURST databases not successfully generated: check input fasta and BURST installation')
	return acx, edx


def align_bgcs(seqs, acx, edx):
	alignment = os.path.join(output, 'bgc_alignments.b6')
	os.system(' '.join(['burst -q', seqs, '-a', acx, '-f -n -r', edx, '-o', alignment]))
	return alignment


def main():
	parser = make_arg_parser()
	args = parser.parse_args()
	outpath = args.output
	bgc_db = args.bgc_fna
	seqs = args.seqs
	if not os.path.isfile(seqs):
		raise ValueError('Sequences file not found')
	threshold = float(args.threshold)
	if not os.path.isdir(outpath):
		os.mkdir(outpath)
		if not os.path.isdir(outpath):
			raise ValueError('Error creating output directory; check given path and try again')
	temp_path = os.path.join(outpath, 'tempz')
	os.mkdir(temp_path)
	acx, edx = create_bgc_db(bgc_db, temp_path)
	alignment = align_bgcs(seqs, acx, edx)

if __name__ == '__main__':
	main()
