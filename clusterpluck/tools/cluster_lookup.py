#!/usr/bin/env Python

import argparse
import sys
import os
import pandas as pd
from collections import defaultdict
from clusterpluck.tools.annotations import refseq_to_name
from clusterpluck.tools.h_clustering import process_hierarchy
from clusterpluck.scripts.run_blastp import run_blastp
from clusterpluck.tools.suppress_print import suppress_stdout
from ninja_dojo.database import RefSeqDatabase
from ninja_dojo.taxonomy import NCBITree
from ninja_utils.parsers import FASTA
from multiprocessing import cpu_count


# The arg parser
def make_arg_parser():
	parser = argparse.ArgumentParser(description='Look up relevant information for given OFU(s), such as DNA and AA sequences, and organism names')
	parser.add_argument('-s', '--scores', help="The appropriate scores matrix resource for this data (csv)", default='-')
	parser.add_argument('-t', '--height', help='The similarity/identity at which the OFUs were picked (0-100)', required=True, type=float)
	parser.add_argument('-m', '--mpfa', help='The .mpfa resource for this ClusterPluck database', required=False)
	parser.add_argument('--mibig', help='To search OFU sequences against the MIBiG database, provide path to the MIBiG blastp database files (including database name, no extension)', required=False)
	parser.add_argument('-d', '--dna_fasta', help='The multi-fasta resource containing cluster DNA sequences for this ClusterPluck database', required=False)
	parser.add_argument('-b', '--bread', help='Where to find the header for the sequence (default="ref|,|")', default='ref|,|')
	parser.add_argument('-o', '--output', help='Directory in which to save the cluster information files (default = cwd)', required=False, default='.')
	parser.add_argument('-c', '--ofu', help='Comma-separated list of the ofus (e.g. ofu00001,ofu00003) on which to provide information', required=False, type=str)
	parser.add_argument('-n', '--name', help='Comma-separated list of the RefSeq IDs (with or without cluster number) for which to provide a list of OFUs', required=False, type=str)
	parser.add_argument('-y', '--types', help='A CSV file containing the predicted product types for each cluster', required=False)
	return parser


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
					ctype = str(ctype.iloc[0, 0])
				else:
					ctype = 'NA'
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


def identify_organism(org, db, nt):
	# with suppress_stdout():
	if 'cluster' in org:
		if '.cluster' in org:
			org.replace('.cluster', '_cluster')
		refseqid = '_'.join(org.split('_')[:2])
	else:
		refseqid = org
	name = refseq_to_name(refseqid, db=db, nt=nt)
	return name


def list_organism_ofus(orgs, hclus, height, outpath):
	bgc_dd = defaultdict(list)
	for value, key in hclus.itertuples(index=True):
		key = str('%05d' % key)
		key = ''.join(['ofu', key])
		bgc_dd[key].extend(list([value]))
	orgs_list = orgs.split(',')
	i = len(orgs_list)
	ofu_dict = defaultdict(list)
	# Preload the Database and Tree
	db = RefSeqDatabase()
	nt = NCBITree()
	for org in orgs_list:
		print(org)
		org_ofu_dup = []
		for ofu_num, ofu_orgs in bgc_dd.items():
			for ofu_org in ofu_orgs:
				if org in ofu_org:
					name = identify_organism(org, db=db, nt=nt)
					if ofu_num not in org_ofu_dup:
						org_ofu_dup.append(ofu_num)
						ofu_dict[name].extend([ofu_num])
					else:
						pass
	height = str(int(height))
	ofu_file = ''.join(['OFUs_from_refseq_id', height, '.txt'])
	outdf = pd.DataFrame.from_dict(ofu_dict, orient='index')
	if not outdf.empty:
		with open(os.path.join(outpath, ofu_file), 'w') as outf:
			outdf.to_csv(outf, sep='\t', header=False)
		print('\nOFU assigments written to file for the %d organism ID entries given.\n' % i)
	else:
		print('\nNo OFU assignments found; check RefSeq ID and format')
	return None


def main():
	parser = make_arg_parser()
	args = parser.parse_args()
	outpath = args.output
	with open(args.scores, 'r') as inf:
		h = 1 - (args.height / 100)
		hclus = process_hierarchy(inf, h)
	if args.types:
		with open(args.types, 'r') as in_t:
			typetable = pd.read_csv(in_t, header=0, index_col=0)
	else:
		typetable = False
	if args.ofu:
		ofus = args.ofu
		bgc_dd = list_organisms(ofus, hclus, typetable, outpath)
		if args.dna_fasta or args.mpfa:  # only generate OFU sequence files if the appropriate files are provided
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
					ofu_aaseqfile = ''.join(['ofu', ofu_n, '_aasequences.mpfa'])
					aa_outf = open(os.path.join(outpath, ofu_aaseqfile), 'w')
					for bgc in bgcs:
						with open(args.mpfa, 'r') as inf_m:
							aa_outf = compile_ofu_sequences(inf_m, bgc, aa_outf)
					aa_outf.close()
					if args.mibig:
						mibigout = ''.join(['ofu', ofu_n, '_vs_MiBiG.txt'])
						mibigout = os.path.join(outpath, mibigout)
						# print(mibigout)
						blastout = open(mibigout, 'w')
						cpus = int(cpu_count() / 2)
						print('Blasting OFU%s against database using %s cpus\n' % (ofu_n, cpus))
						mibig_db = str(os.path.relpath(args.mibig))
						# if ' ' in mibig_db:
						# 	mibig_db = mibig_db.replace(' ', '\ ')
						# print(mibig_db)
						ofu_query = str(os.path.join(outpath, ofu_aaseqfile))
						# if ' ' in ofu_query:
						# 	ofu_query = ofu_query.replace(' ', '\ ')
						# print(ofu_query)
						blastresult = run_blastp(ofu_query, mibigout, mibig_db, cpus)
						blastout.write(blastresult)
						blastout.close()
				if args.dna_fasta:
					ofu_dnaseqfile = ''.join(['ofu', ofu_n, '_dnasequences.fna'])
					dna_outf = open(os.path.join(outpath, ofu_dnaseqfile), 'w')
					for bgc in bgcs:
						with open(args.dna_fasta, 'r') as inf_d:
							dna_outf = compile_ofu_dnasequences(inf_d, bgc, dna_outf)
					dna_outf.close()
			print('Sequence files written for %d OFUs.\n' % i)
		else:
			pass
	if args.name:  # if using to look up OFUs from the RefSeq IDs, run this bit
		orgs = args.name
		height = args.height
		list_organism_ofus(orgs, hclus, height, outpath)

	sys.exit()


if __name__ == '__main__':
	main()
