#!/usr/bin/env python

# Run from antismash_results directory!!
# Parses out the predicted cluster types from the antismash output file.
# With -c flag, also compiles all the cluster types for all genome results present
# Returns CSV files with the refseq accession and the cluster types

import argparse
import os
import os.path
import re
import pandas as pd
import csv

from collections import defaultdict
from ninja_dojo.database import RefSeqDatabase
from ninja_dojo.taxonomy import NCBITree
from clusterpluck.tools.annotations import refseq_to_name
from clusterpluck.tools.annotations import refseq_to_tid

usage = 'list_cluster_types.py'


def make_arg_parser():
	parser = argparse.ArgumentParser(description='Lists the predicted product type of each cluster in a csv')
	parser.add_argument('-c', '--compile', help='Also compile all cluster types for each genome.', action='store_true', default=False)
	parser.add_argument('-a', '--annotate', help='After compiling, annotate the index with NCBI TID, RefSeq Accession, and Organism name', action='store_true', default=False)
	return parser


def parse_products(rdir):
	filelist = os.listdir(rdir)
	# debug
	# print(filelist)
	# define the start of the sequence by the CDS line
	title_begin = False
	sequence_begin = False
	if 'cluster_types' not in filelist:
		os.mkdir(os.path.join(rdir, 'cluster_types'))
	i = 0
	dd = defaultdict(dict)
	for infilename in filelist:
		if infilename.endswith('.gbk'):
			if infilename.endswith('final.gbk'):
				pass
			# elif not infilename.startswith('NC'):
			# 	pass
			else:
				# print(infilename)
				i += 1
				clusterid = infilename.replace('.gbk', '')
				clusterid = clusterid.replace('.clu', '_clu')
				with open(os.path.join(rdir, infilename), 'r') as infile:
					for line in infile:
						if line.startswith("                     /product"):
							p = re.compile(r"^(\s+)\/(product)=\"(.*)\"$")
							m = p.search(line)  # searches using the regex defined above
							prod = m.group(3)
							dd[clusterid] = prod
						else:
							pass
				infile.close()
	if i == 0:
		os.rmdir(os.path.join(rdir, 'cluster_types'))
	else:
		outfname = '_'.join(clusterid.split('_')[0:2])
		outfile = open(os.path.join(rdir, 'cluster_types', outfname + '.csv'), 'w')
		cdf = pd.DataFrame.from_dict(dd, orient='index')
		cdf.columns = ['cluster_type']
		cdf.sort_index(axis=0, inplace=True)
		outfile.write(cdf.to_csv())
		outfile.close()


def compile_types(cdir, outf):
	for file in os.listdir(os.path.join(cdir, 'cluster_types')):
		if file.endswith('.csv'):
			# print(file)
			with open(os.path.join(cdir, 'cluster_types', file), 'r') as infile:
				next(infile)
				for line in infile:
					outf.write(line)
			infile.close()
		else:
			pass
	return outf


def main():
	parser = make_arg_parser()
	args = parser.parse_args()
	if os.getcwd().split('/')[-1] == 'antismash_results':
		for rdir in os.listdir('.'):
			if rdir.startswith('GCF'):
				parse_products(rdir)
			else:
				pass
	else:
		print('\nERROR:\nYou must run this script from within the antismash_results directory containing the results for each genome\n')
		quit()
	if args.compile:
		if "compiled_cluster_types" not in os.listdir('.'):
			os.mkdir("compiled_cluster_types")
		outfilename = 'compiled_cluster_types.csv'
		with open(os.path.join("compiled_cluster_types", outfilename), 'w') as outf:
			for cdir in os.listdir('.'):
				if cdir.startswith('GCF'):
					if "cluster_types" not in os.listdir(cdir):
						pass
					else:
						outf = compile_types(cdir, outf)
			outf.close()
	if args.annotate:
		if not args.compile:
			print('\nSORRY:\nAnnotation only available with compiled option (-c)\n')
			quit()
		else:
			# Preload the Database and Tree
			db = RefSeqDatabase()
			nt = NCBITree()
			strain_label = []
			with open(os.path.join('compiled_cluster_types', 'compiled_cluster_types.csv')) as intab:
				odf = pd.read_csv(intab, index_col=0)
				refseq_list = list(odf.index)
				for refseq_id in refseq_list:
					organism = refseq_to_name(refseq_id, db=db, nt=nt)
					ncbi_tid = refseq_to_tid(refseq_id, db=db)
					ncbi_tid = str(ncbi_tid)
					# genus_species = organism.split(';')[-1]
					# genus_species = genus_species.replace('s__', '')
					if ncbi_tid == organism:  # sometimes DOJO can't look up the refseq accession; in this case, just return refseq.
						strain_label.append(refseq_id)
					else:
						strain_label.append('ncbi_tid|%s|ref|%s|organism|%s' % (ncbi_tid, refseq_id, organism))
				odf.index = strain_label
				an_outn = 'annotated_cluster_types.csv'
			with open(os.path.join('compiled_cluster_types', an_outn), 'w') as an_outf:
				odf.to_csv(an_outf)
	else:
		pass


if __name__ == '__main__':
	main()
