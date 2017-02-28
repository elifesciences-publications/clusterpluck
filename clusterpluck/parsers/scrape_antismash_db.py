#!/usr/bin/env python

# Run from directory above all the database gbks

import argparse
import os
import re
import csv
import sys
import pandas as pd
import logging
from collections import defaultdict
from dojo.taxonomy import NCBITree
from ninja_utils.parsers import FASTA


# The arg parser
def make_arg_parser():
	parser = argparse.ArgumentParser(description='Description here')
	parser.add_argument('-i', '--input', help='Directory in which to find all the cluster .gbk files (default = cwd)', required=False, default='.')
	parser.add_argument('-nt_cat', help='The nucleotide_catalog tsv file matching Genbank accessions to NCBI taxon IDs', required=True)
	parser.add_argument('-o', '--output', help='Where to save the output files (default = new dir in cwd)', required=False, default='.')
	parser.add_argument('--no_compile', help='Do not compile all the cluster files and DNA files into single .mpfa and .fna files',
						action='store_true', required=False, default=False)
	parser.add_argument('--just_compile', help='Only run the compilation step - specify the output locations with -o flag.',
						action='store_true', required=False, default=False)
	return parser


def parse_aa_seqs(gbk_file, tid_org, gbk_filepath, outpath):
	ncbi_tid = str(tid_org[0])
	organism = str(tid_org[1])
	header = '>ncbi_tid|%s|genbank|' % ncbi_tid  # start each sequence's identifier with the pacman >
	# define the start of the sequence by the CDS line
	title_begin = False
	sequence_begin = False
	mpfa_results = 'antismash_db_protein_seqs'
	if mpfa_results not in os.listdir(outpath):
		os.mkdir(os.path.join(outpath, mpfa_results))
	i = 0
	if gbk_file.endswith('.gbk'):
		if gbk_file.endswith('final.gbk'):
			pass
		else:
			# print(infilename)
			i += 1
			header_c = header + gbk_file.replace('.gbk', '')
			header_c = header_c.replace('.clu', '_clu')
			outfilename = 'aa_' + gbk_file + '.mpfa'
			outfilename = outfilename.replace('.gbk', '').replace('.clu', '_clu')
			outfile = open(os.path.join(outpath, mpfa_results, outfilename), 'w')
			with open(os.path.join(gbk_filepath), 'r') as infile:
				for line in infile:
					if title_begin:  # only do this if 'CDS  ' starts the line
						if line.startswith("                     /locus_tag"):
							p = re.compile(r"^(\s+)(\/locus_tag=)\"(.*)\"")
							m = p.search(line)  # searches using the regex defined above
							outfile_m = m.group(3)
							outfile.write(header_c + '_' + outfile_m)  # use the filename to ID the file on the first line
							outfile.write('|organism|%s\n' % organism)
						if line.startswith('                     /translation'):
							sequence_begin = True
						if sequence_begin:
							if line.startswith('                     /translation'):
								aa_p = re.compile(r"^(\s+)(\/translation=\")([A-Z]+)")
								aa_m = aa_p.search(line)  # searches using the regex defined above
								out_aa = aa_m.group(3)
								outfile.write(out_aa)
							if line.startswith('                     '):
								outfile.write(''.join([ch for ch in line if ch in set('G,A,L,M,F,W,K,Q,E,S,P,V,I,C,Y,H,R,N,D,T')]))
							else:
								outfile.write('\n')
								sequence_begin = False
								if line.startswith('     CDS  '):
									title_begin = True
								else:
									title_begin = False
									sequence_begin = False
					elif line.startswith('     CDS  '):
						title_begin = True  # identifies the line starting with CDS as cluster module sequence start
			outfile.close()
	return None


def parse_dna_seqs(gbk_file, tid_org, gbk_filepath, outpath):
	ncbi_tid = str(tid_org[0])
	organism = str(tid_org[1])
	header = '>ncbi_tid|%s|genbank|' % ncbi_tid  # start each sequence's identifier with the pacman >
	# define the start of the sequence by the ORIGIN line
	sequence_begin = False
	dna_results = 'antismash_db_dna_seqs'
	if dna_results not in os.listdir(outpath):
		os.mkdir(os.path.join(outpath, dna_results))
	i = 0
	if gbk_file.endswith('.gbk'):
		if gbk_file.endswith('final.gbk'):
			pass
		else:
			# print(infilename)
			i += 1
			header_c = header + gbk_file.replace('.gbk', '')
			header_c = header_c.replace('.clu', '_clu')
			outfilename = 'dna_' + gbk_file + '.fna'
			outfilename = outfilename.replace('.gbk', '').replace('.clu', '_clu')
			outfile = open(os.path.join(outpath, dna_results, outfilename), 'w')
			outfile.write(header_c + '|organism|%s\n' % organism)
			with open(os.path.join(gbk_filepath), 'r') as infile:
				for line in infile:
					if sequence_begin:  # only do this if ORIGIN starts the line
						# joins together only the characters on the line in the set atcg
						outfile.write(''.join([ch for ch in line if ch in set(('a', 't', 'c', 'g'))]))
					elif line.startswith('ORIGIN'):
						sequence_begin = True  # identifies the line starting with ORIGIN as sequence start
			outfile.close()
	return None


def parse_cluster_types(gbkpath, outpath, gbk_dd):
	filelist = os.listdir(gbkpath)
	# debug
	# print(filelist)
	# define the start of the sequence by the CDS line
	cluster_begin = False
	if 'antismash_db_product_types' not in os.listdir(outpath):
		os.mkdir(os.path.join(outpath, 'antismash_db_product_types'))
	i = 0
	type_dd = defaultdict(dict)
	for gbk in filelist:
		if gbk.endswith('.gbk'):
			if gbk.endswith('final.gbk'):
				pass
			# print(infilename)
			i += 1
			clusterid = gbk.replace('.gbk', '')
			clusterid = clusterid.replace('.clu', '_clu')
			gbk_id = gbk.split('.cluster')[0]
			# tid_org = []
			tid_org = gbk_dd[gbk_id]
			if not tid_org:
				tid_org = ['na', 'k__None;p__None;c__None;o__None;f__None;g__None;s__None;t__None']
			ncbi_tid = str(tid_org[0])
			organism = str(tid_org[1])
			cluster_label = 'ncbi_tid|%s|genbank|%s|organism|%s' % (ncbi_tid, clusterid, organism)
			with open(os.path.join(gbkpath, gbk), 'r') as in_gbk:
				for line in in_gbk:
					if cluster_begin:
						if line.startswith("                     /product"):
							if line.endswith("\"\n"):
								p = re.compile(r"^(\s+)\/(product)=\"(.*)\"")
								m = p.search(line)
								prod = str(m.group(3))
							else:
								p = re.compile(r"^(\s+)\/(product)=\"(.*)$")
								m = p.search(line)
								nxline = next(in_gbk)
								p2 = re.compile(r"(                     )(.*)\"")
								m2 = p2.search(nxline)
								prod = ' '.join([str(m.group(3)), str(m2.group(2))])
							type_dd[cluster_label] = prod
						elif line.startswith("     gene"):
							cluster_begin = False
					elif line.startswith("     cluster"):
						cluster_begin = True
	cdf = pd.DataFrame.from_dict(type_dd, orient='index')
	cdf.columns = ['cluster_type']
	cdf.sort_index(axis=0, inplace=True)
	with open(os.path.join(outpath, 'antismash_db_product_types', 'asdb_product_types.csv'), 'w') as outfile:
		cdf.to_csv(outfile)
	return None


def compile_files(outpath):
	protein_seqs = os.path.join(outpath, 'asDB_protein_seqs.mpfa')
	dna_seqs = os.path.join(outpath, 'asDB_dna_seqs.fna')
	with open(protein_seqs, 'w') as aa_outfile:
		for aafile in os.listdir(os.path.join(outpath, 'antismash_db_protein_seqs')):
			aafile = os.path.join(outpath, 'antismash_db_protein_seqs', aafile)
			with open(aafile, 'r') as aa_in:
				fasta_gen = FASTA(aa_in)
				for header, sequence in fasta_gen.read():
					aa_outfile.write('>' + header + '\n')
					aa_outfile.write('>' + sequence + '\n')
	aa_outfile.close()
	with open(dna_seqs, 'w') as dna_outfile:
		for dnafile in os.listdir(os.path.join(outpath, 'antismash_db_dna_seqs')):
			dnafile = os.path.join(outpath, 'antismash_db_dna_seqs', dnafile)
			with open(dnafile, 'r') as dna_in:
				fasta_gen = FASTA(dna_in)
				for header, sequence in fasta_gen.read():
					dna_outfile.write('>' + header + '\n')
					dna_outfile.write('>' + sequence + '\n')
	dna_outfile.close()
	return None


def tid_to_name(tid, nt=NCBITree()):
	tid = int(tid)
	organism = nt.green_genes_lineage(tid, depth=8, depth_force=True)
	return organism


def main():
	parser = make_arg_parser()
	args = parser.parse_args()
	nt_cat = os.path.join(args.nt_cat)
	gbkpath = os.path.join(args.input)
	outpath = os.path.join(args.output)
	if args.just_compile:
		compile_files(outpath)
		sys.exit()
	if not os.path.isdir(outpath):
		os.mkdir(os.path.join(outpath))
		if not os.path.isdir(outpath):
			print('\nError creating output directory; check given path and try again\n')
			sys.exit()
	logfile = os.path.join(outpath, 'scrapelog.log')
	logging.basicConfig(filename=logfile, level=logging.DEBUG, format='%(asctime)s %(message)s', datefmt='%m/%d/%Y %I:%M:%S %p')
	gbks = os.listdir(gbkpath)
	gbks = [f for f in gbks if f.endswith('gbk')]
	with open(nt_cat, 'r') as nt_catalog:
		gbk_dd = defaultdict(list)
		reader = csv.reader(nt_catalog, delimiter='\t')
		next(reader)
		nt = NCBITree()
		gbk_set = set()
		for gbk_file in gbks:
			gbk_id = gbk_file.split('.cluster')[0]
			gbk_set.add(gbk_id)
		for line in reader:
			if line[1] in gbk_set:
				tid = line[2]
				organism = tid_to_name(tid, nt=nt)
				# print(line[1] + tid + organism)
				gbk_dd[line[1]] = [tid, organism]
	i = 0
	for gbk_file in gbks:
		gbk_id = gbk_file.split('.cluster')[0]
		tid_org = gbk_dd[gbk_id]
		if not tid_org:
			print('Error getting taxonomy for %s for cluster file %s' % (gbk_id, gbk_file))
			logging.warning('Error getting taxonomy for %s for cluster file %s' % (gbk_id, gbk_file))
			tid_org = ['na', 'k__None;p__None;c__None;o__None;f__None;g__None;s__None;t__None']
			i += 1
		# print(tid_org)
		# ncbi_tid = str(tid_org[0])
		# organism = str(tid_org[1])
		gbk_filepath = os.path.join(gbkpath, gbk_file)
		parse_aa_seqs(gbk_file, tid_org, gbk_filepath, outpath)
		parse_dna_seqs(gbk_file, tid_org, gbk_filepath, outpath)
	parse_cluster_types(gbkpath, outpath, gbk_dd)
	if not args.no_compile:
		compile_files(outpath)
	logging.warning('DOJO could not acquire NCBI tid information for %s clusters' % i)

if __name__ == '__main__':
	main()
