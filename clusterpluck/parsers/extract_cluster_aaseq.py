#!/usr/bin/env python

# Run from antismash_results directory!!
# Parses .gbk output files of antiSMASH 3.0, yielding a txt file with each module's
# coding sequence (CDS) listed in an faa-like file. The cluster module ID is listed as:
# >NC_011593.1.cluster001_ctg1_orf00220 for example.

import os
import os.path
import re
from clusterpluck.tools.compile_mpfa import compile_aa

usage = 'extract_cluster_aaseq.py'

def make_arg_parser():
	parser = argparse.ArgumentParser(description='Convert blastp output txt table to a scores matrix in csv format')
	parser.add_argument('-c', '--compile', help='Also compile all ORF amino acid sequences per genome.', action='store_true', default=False)
	return parser

def parse_aa(rdir):
	filelist = os.listdir(rdir)
	# debug
	# print(filelist)
	header = '>'  # start each sequence's identifier with the pacman >
	# define the start of the sequence by the CDS line
	codingstart = 'CDS   '
	title_begin = False
	sequence_begin = False
	if 'cluster_aa_sequences' not in filelist:
		os.mkdir(os.path.join(rdir, 'cluster_aa_sequences'))
	i = 0
	for infilename in filelist:
		if infilename.endswith('.gbk'):
			if infilename.endswith('final.gbk'):
				pass
			else:
				i += 1
				header_f = header + infilename.replace('.gbk', '')
				header_f = header_f.replace('.clu', '_clu')
				outfilename = 'aa_' + infilename + '.mpfa'
				outfilename = outfilename.replace('.gbk', '')
				outfile = open(os.path.join(rdir, 'cluster_aa_sequences', outfilename), 'w')
				with open(os.path.join(rdir, infilename), 'r') as infile:
					for line in infile:
						if title_begin:  # only do this if 'CDS  ' starts the line
							if line.startswith("                     /locus_tag"):
								p = re.compile(r"^(\s+)(/locus_tag=)\"(ctg)(\d_\w+)\"")
								m = p.search(line)  # searches using the regex defined above
								outfileM = ''.join(m.group(3, 4))
								outfile.write(header_f + '_')  # use the filename to ID the file on the first line
								outfile.write(outfileM + '\n')
							if line.startswith('                     /translation'):
								sequence_begin = True
							if sequence_begin:
								if line.startswith('                     /translation'):
									aa_p = re.compile(r"^(\s+)(\/translation\=\")([A-Z]+)")
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
			# 							if line.startswith('                     /translation'):
			# 								sequence_begin = True
									else:
										title_begin = False
										sequence_begin = False
						elif line.startswith('     CDS  '):
							title_begin = True  # identifies the line starting with CDS as cluster module sequence start
		# else:
	if i == 0:
		os.rmdir(os.path.join(rdir, 'cluster_aa_sequences'))
	else:
		infile.close()
		outfile.close()


def main():
	parser = make_arg_parser()
	args = parser.parse_args()
	if os.getcwd().split('/')[-1] == 'antismash_results':
		for rdir in os.listdir('.'):
			if rdir.startswith('GCF'):
				parse_aa(rdir)
			else:
				pass
	else:
		print('\nERROR:\nYou must run this script from within the antismash_results directory containing the results for each genome\n')
		quit()

	if args.compile:
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
