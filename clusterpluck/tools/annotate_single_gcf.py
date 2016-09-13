#!/usr/bin/env python
import argparse
import sys
import os

from ninja_utils.parsers import FASTA
from ninja_dojo.database import RefSeqDatabase


# The arg parser for this wrapper
def make_arg_parser():
    parser = argparse.ArgumentParser(description='Reformat fasta header to include refseq, assembly, and ncbi_tid')
    parser.add_argument('-i', '--input', help='Input is a FASTA file.', default='-')
    parser.add_argument('-o', '--output', help='If nothing is given, then stdout, else write to file', default='-')
    parser.add_argument('-v', '--verbose', help='Print extra statistics', action='store_true', default=False)
    return parser



def main():
    parser = make_arg_parser()
    args = parser.parse_args()

    db = RefSeqDatabase()
    # parse command line
    with open(args.input, 'r') if args.input != '-' else sys.stdin as inf:
        fasta_gen = FASTA(inf)
        assembly_version = os.path.basename(args.input).split('_genomic')[0]
        with open(args.output, 'w') if args.output != '-' else sys.stdout as outf:
            for header, sequence in fasta_gen.read():
                ncbi_tid = db.get_ncbi_tid_from_refseq_accession_version(header.split()[0])[0]
                outf.write('>ncbi_tid|%d|assembly|%s|ref|%s\n' % (ncbi_tid, assembly_version, header))
                outf.write(sequence+'\n')

if __name__ == '__main__':
    main()
