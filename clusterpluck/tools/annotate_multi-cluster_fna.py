#!/usr/bin/env python
import argparse
import sys
import os

from ninja_utils.parsers import FASTA
from ninja_dojo.database import RefSeqDatabase
from ninja_dojo.taxonomy import NCBITree

# The arg parser for this wrapper
def make_arg_parser():
	parser = argparse.ArgumentParser(description='Reformat multi-fasta file header to include refseq id, ncbi_tid, and Genus_species')
	parser.add_argument('-i', '--input', help='Input is a multi-cluster FASTA file.', default='-')
	parser.add_argument('-o', '--output', help='If nothing is given, then stdout, else write to file', default='-')
	parser.add_argument('-v', '--verbose', help='Print extra statistics', action='store_true', default=False)
	return parser


def main():
	parser = make_arg_parser()
	args = parser.parse_args()

	db = RefSeqDatabase()
	nt = NCBITree()
	# parse command line
	with open(args.input, 'r') if args.input != '-' else sys.stdin as inf:
		fasta_gen = FASTA(inf)
		assembly_version = os.path.basename(args.input).split('_genomic')[0]
		with open(args.output, 'w') if args.output != '-' else sys.stdout as outf:
			for header, sequence in fasta_gen.read():
				if '.cluster' in header:
					header = header.replace('.cluster','_cluster')
				else:
					pass
				ncbi_tid = db.get_ncbi_tid_from_refseq_accession_version(header.split('_cluster')[0])
				if ncbi_tid:
					ncbi_tid = ncbi_tid[0]
					organism = nt.gg_lineage(ncbi_tid)
					genus_species = organism.split(';')[-1]
					genus_species = genus_species.replace('s__','')
					outf.write('>ncbi_tid|%d|ref|%s|organism|%s|\n' % (ncbi_tid, header, genus_species))
					outf.write(sequence+'\n')
				else:
					outf.write('>ref|%s|\n' % (header))
					outf.write(sequence+'\n')

if __name__ == '__main__':
	main()
