#!/usr/bin/env Python

import argparse
import sys
import os
import pandas as pd
from collections import defaultdict
from clusterpluck.tools.annotations import refseq_to_name
from clusterpluck.tools.h_clustering import process_hierarchy
from ninja_dojo.database import RefSeqDatabase
from ninja_dojo.taxonomy import NCBITree
from ninja_utils.parsers import FASTA
from contextlib import contextmanager

# Take a list of strain_clusters (in the format "NC_004307.2_cluster004") and search the .mpfa database file
# to return a protein fasta (.mpfa) file containing the sequences of the clusters requested.
#
# >ncbi_tid|206672|ref|NC_004307.2_cluster004_ctg1_orf02030|organism|Bifidobacterium_longum|


# The arg parser
def make_arg_parser():
	parser = argparse.ArgumentParser(description='Look up relevant information for given OFU(s), such as DNA and AA sequences, and organism names')
	parser.add_argument('-s', '--scores', help="The appropriate scores matrix resource for this data (csv)", default='-')
	parser.add_argument('-t', '--height', help='The similarity/identity at which the OFUs were picked (0-100)', required=True, type=float)
	parser.add_argument('-m', '--mpfa', help='The .mpfa resource for this ClusterPluck database', required=False)
	parser.add_argument('-d', '--dna_fasta', help='The multi-fasta resource containing cluster DNA sequences for this ClusterPluck database', required=False)
	parser.add_argument('-b', '--bread', help='Where to find the header for the sequence (default="ref|,|")', default='ref|,|')
	parser.add_argument('-o', '--output', help='Directory in which to save the cluster information files (default = cwd)', required=False, default='.')
	parser.add_argument('-c', '--ofu', help='Comma-separated list of the ofus (e.g. ofu00001,ofu00003) on which to provide information', required=False, type=str)
	parser.add_argument('-y', '--types', help='A CSV file containing the predicted product types for each cluster', required=False)
	return parser


@contextmanager
def suppress_stdout():
	with open(os.devnull, 'w') as devnull:
		old_stdout = sys.stdout
		sys.stdout = devnull
		try:
			yield
		finally:
			sys.stdout = old_stdout


def list_organisms(ofus, hclus, typetable, outpath):
	bgc_dd = defaultdict(list)
	for value, key in hclus.itertuples(index=True):
		key = str('%05d' % key)
		bgc_dd[key].extend(list([value]))
	ofu_list = ofus.split(',')
	i = 0
	# Preload the Database and Tree
	db = RefSeqDatabase()
	nt = NCBITree()
	for ofu in ofu_list:
		ofu = str(ofu)
		if ofu.startswith('ofu'):
			ofu_n = str(ofu.replace('ofu', ''))
		else:
			ofu_n = ofu
		bgcs = bgc_dd[ofu_n]
		name_dict = defaultdict(list)
		with suppress_stdout():
			for bgc in bgcs:
				refseqid = '_'.join(bgc.split('_')[:2])
				name = refseq_to_name(refseqid, db=db, nt=nt)
				if typetable is not False:
					ctype = typetable.filter(like=bgc, axis=0)
					ctype = str(ctype.iloc[0,0])
				if bgc == name:
					name_dict[bgc] = [ctype, refseqid]
				else:
					name_dict[bgc] = [ctype, name]
		ofu_file = ''.join(['ofu', ofu_n, '.txt'])
		with open(os.path.join(outpath, ofu_file), 'w') as outf:
			outdf = pd.DataFrame.from_dict(name_dict, orient='index')
			outdf.columns = ['predicted_type', 'organism']
			outdf.to_csv(outf, sep='\t')
		i += 1
	print('\nOrganism information for %d OFUs written to file.\n' % i)
	return bgc_dd


def compile_ofu_sequences(inf_m, bgc, aa_outf):
	mpfa_gen = FASTA(inf_m)
	bgc = str(bgc)
	for header, sequence in mpfa_gen.read():
		if '.cluster' in header:
			header = header.replace('.cluster', '_cluster')
		if bgc in header:
			aa_outf.write(''.join(['>', header, '\n']))
			aa_outf.write(''.join([sequence, '\n']))
	return aa_outf


def compile_ofu_dnasequences(inf_d, bgc, dna_outf):
	fasta_gen = FASTA(inf_d)
	bgc = str(bgc)
	for header, sequence in fasta_gen.read():
		if '.cluster' in header:
			header = header.replace('.cluster', '_cluster')
		if bgc in header:
			dna_outf.write(''.join(['>', header, '\n']))
			dna_outf.write(''.join([sequence, '\n']))
	return dna_outf


def main():
	parser = make_arg_parser()
	args = parser.parse_args()
	outpath = args.output
	with open(args.scores, 'r') as inf:
		h = 1 - (args.height / 100)
		hclus = process_hierarchy(inf, h)
	ofus = args.ofu
	if args.types:
		with open(args.types, 'r') as in_t:
			typetable = pd.read_csv(in_t, header=0, index_col=0)
	else:
		typetable = False
	bgc_dd = list_organisms(ofus, hclus, typetable, outpath)
	if args.dna_fasta or args.mpfa:
		ofu_list = ofus.split(',')
		i = 0
		for ofu in ofu_list:
			i += 1
			if ofu.startswith('ofu'):
				ofu_n = str(ofu.replace('ofu', ''))
			else:
				ofu_n = ofu
			bgcs = bgc_dd[ofu_n]
			if args.mpfa:
				ofu_aaseqfile = ''.join(['ofu', ofu_n, '_aasequences.txt'])
				aa_outf = open(os.path.join(outpath, ofu_aaseqfile), 'w')
				for bgc in bgcs:
					with open(args.mpfa, 'r') as inf_m:
						aa_outf = compile_ofu_sequences(inf_m, bgc, aa_outf)
				aa_outf.close()
			if args.dna_fasta:
				ofu_dnaseqfile = ''.join(['ofu', ofu_n, '_dnasequences.txt'])
				dna_outf = open(os.path.join(outpath, ofu_dnaseqfile), 'w')
				for bgc in bgcs:
					with open(args.dna_fasta, 'r') as inf_d:
						dna_outf = compile_ofu_dnasequences(inf_d, bgc, dna_outf)
				dna_outf.close()
		print('Sequence files written for %d OFUs.\n' % i)
	else:
		pass
	sys.exit()


if __name__ == '__main__':
	main()
