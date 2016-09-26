#!/usr/bin/env Python

import argparse
import sys
import csv
import numpy as np
import pandas as pd
from collections import defaultdict
from clusterpluck.tools.annotations import refseq_to_tid
from clusterpluck.tools.annotations import refseq_to_name

from ninja_dojo.database import RefSeqDatabase
from ninja_dojo.taxonomy import NCBITree


# The arg parser
def make_arg_parser():
	parser = argparse.ArgumentParser(description='Generates a binary OFU matrix (genomes vs BGCs) from clustered clusters')
	parser.add_argument('-i', '--input', help='Input file: a CSV of the BGCs and hclust clusters from R.', required=True)
	parser.add_argument('-o', '--output', help='Where to save the output csv; default to screen', required=False, default='-')
	parser.add_argument('-a', '--annotate', help='Annotate the OFU table with NCBI tid, RefSeq Accession, and organism name', action='store_true', default=False)
	return parser


# initialize the dict with zeros for each OFU
def outer(size):
	return lambda: np.zeros(size, dtype=np.int)


# build the ofu table from an hclus (in R) generated csv
def cluster_ofus(inf2, dd):
	ofu_names = []
	hcsv = csv.reader(inf2, delimiter=',')
	next(hcsv, None)
	for line in hcsv:
		name = line[0]
		refseq_id = '_'.join(name.split('_')[0:2])
		ofu = int(line[1])
		ofu_name = ('ofu', str('%03d' % int(ofu)))
		ofu_names.append('_'.join(ofu_name))
		dd[refseq_id][ofu] += 1  # adds a "1" to the reference for that clustered OFU for this organism
	df = pd.DataFrame.from_dict(dd)
	df = df.T
	df.drop([0], axis=1, inplace=True)  # removes the empty first ("0"th) OFU column
	ofu_cols = list(set(ofu_names))
	ofu_cols.sort()
	df.columns = ofu_cols
	df.sort_index(axis=0)
	df.sort_index(axis=1)
	return df


def main():
	parser = make_arg_parser()
	args = parser.parse_args()
	# Parse command line

	with open(args.input, 'r') as inf:
		hclus = pd.read_csv(inf, sep=',', header=0, index_col=0)
		size = hclus.max(0)[0]  # get the total number of clustered OFUs at the height cutoff used in R
		size += 1
		fill = outer(size)
		dd = defaultdict(fill)
		with open(args.input, 'r') as inf2:
			df = cluster_ofus(inf2, dd)
			if args.annotate:
				# Preload the Database and Tree
				db = RefSeqDatabase()
				nt = NCBITree()
				strain_label = []
				refseq_list = list(df.index)
				for refseq_id in refseq_list:
					organism = refseq_to_name(refseq_id, db=db, nt=nt)
					ncbi_tid = refseq_to_tid(refseq_id, db=db)
					ncbi_tid = str(ncbi_tid)
					genus_species = organism.split(';')[-1]
					genus_species = genus_species.replace('s__', '')
					strain_label.append('ncbi_tid|%s|ref|%s|organism|%s' % (ncbi_tid, refseq_id, genus_species))
				df.index = strain_label
			else:
				pass
		with open(args.output, 'w') if args.output != '-' else sys.stdout as outf:
			df.to_csv(outf)


if __name__ == '__main__':
	main()
