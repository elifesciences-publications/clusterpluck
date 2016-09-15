#!/usr/bin/env Python

import os
import os.path

usage = "execute from antismash_results directory to compile each genome's cluster into one .faa file"


def compile_aa(cdir, outfile):
	for file in os.listdir(os.path.join(cdir, 'cluster_aa_sequences/')):
		if file.endswith('.mpfa'):
			with open(os.path.join(cdir, 'cluster_aa_sequences/', file), 'r') as infile:
				for line in infile:
					outfile.write(line)
	return outfile


def main():
	if "compiled_cluster_aa_seqs" not in os.listdir('.'):
		os.mkdir("compiled_cluster_aa_seqs")
	for cdir in os.listdir('.'):
		if cdir.startswith('GCF'):
			if "cluster_aa_sequences" not in os.listdir(cdir):
				pass
			else:
				fname = cdir.strip('_genomic')
				outfilename = fname + '_cluster_aa_seqs.mpfa'
				outfile = open(os.path.join("compiled_cluster_aa_seqs", outfilename), 'w')
				compile_aa(cdir, outfile)
				outfile.close()
		else:
			pass


if __name__ == '__main__':
	main()
