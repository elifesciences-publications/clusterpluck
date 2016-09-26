#!/usr/bin/env python

import argparse
import sys

from ninja_dojo.database import RefSeqDatabase
from ninja_dojo.taxonomy import NCBITree

usage = 'gcf_demystifier.py -i GCF_XXXXXXXXXX.X [-o] OUTPUT_FILE [-v]'


# The arg parser for this wrapper
def make_arg_parser():
	parser = argparse.ArgumentParser(description='Annotate from refseq id to tid, refseq, and/or organism taxonomy')
	parser.add_argument('-a', '--assembly', help='Input is a GCF', default='-')
	parser.add_argument('-r', '--refseq', help='Input is the RefSeq Accession ID', default='-')
	parser.add_argument('-t', '--tid', help='Input is the ncbi tid', default='-')
	parser.add_argument('-o', '--output', help='If nothing is given, then stdout, else write to file', default='-')
	parser.add_argument('-v', '--verbose', help='Print extra statistics', action='store_true', default=False)
	return parser


def assembly_to_tid(assembly):
	db = RefSeqDatabase()
	ncbi_tid = db.get_ncbi_tid_from_assembly_accession_version(assembly)[0]
	return ncbi_tid


def refseq_to_tid(refseq_id, db = RefSeqDatabase()):
	ncbi_tid = db.get_ncbi_tid_from_refseq_accession(refseq_id)[0]
	return ncbi_tid


def refseq_to_name(refseq_id, db = RefSeqDatabase(), nt = NCBITree()):
	ncbi_tid = db.get_ncbi_tid_from_refseq_accession(refseq_id)[0]
	organism = nt.green_genes_lineage(ncbi_tid)
	return organism


def main():
	parser = make_arg_parser()
	args = parser.parse_args()

	# parse command line
	with open(args.output, 'w') if args.output != '-' else sys.stdout as outf:
		if args.assembly != '-':
			ncbi_tid = assembly_to_tid(assembly=args.assembly)
		elif args.refseq != '-':
			ncbi_tid = refseq_to_tid(refseq_id=args.refseq)
			organism = refseq_to_name(refseq_id=args.refseq)
			genus_species = organism.split(';')[-1]
			genus_species = genus_species.replace('s__','')
		elif args.tid != '-':
			ncbi_tid = int(args.tid)
		outf.write('>ncbi_tid|%d|organism|%s\n' % (ncbi_tid, genus_species))
		outf.write('\n')

if __name__ == '__main__':
	main()

