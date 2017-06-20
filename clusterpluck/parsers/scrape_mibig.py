#!/usr/bin/env python

# Run from directory above all the database gbks

import argparse
import os
import re
import sys
from collections import defaultdict
from dojo.taxonomy import NCBITree
from ninja_utils.parsers import FASTA


# The arg parser
def make_arg_parser():
	parser = argparse.ArgumentParser(description='Modifies the MIBiG protein fasta file to play well with clusterpluck')
	parser.add_argument('-i', '--input', help='Directory in which to find all the MIBiG .gbk files (default = cwd)', required=True, default='.')
	parser.add_argument('--mibig_aa', help='The existing MIBiG protein fasta, to modify', required=False)
	parser.add_argument('-o', '--output', help='Where to save the output files (default = new dir in cwd)', required=False, default='.')
	return parser


def get_tids(gbk, gbk_file, outf, nt=NCBITree()):
	dict_list = []
	with open(gbk_file, 'r') as inf:
		p = re.compile(r"^(\s+)\/db_xref=\"taxon:(\d+)\"")
		for line in inf:
			if line.startswith('                     /db_xref="taxon:'):
				m = p.search(line)
				ncbi_tid = int(m.group(2))
				organism = nt.green_genes_lineage(ncbi_tid, depth=8, depth_force=True)
				dict_list = [ncbi_tid, organism]
				outf.write(gbk + '\t' + str(ncbi_tid) + '\t' + organism + '\n')
				break
		if not dict_list:
			print(gbk + ' failed to find tid')
			return ['None', 'None']
		else:
			return dict_list


def pluckify_mibig(inf, outaa, gbk_dd):
	mibig_orig = FASTA(inf)
	for header, sequence in mibig_orig.read():
		m_head = header.split('|')
		# print(m_head)
		bgc_id = m_head[0]
		bgc_info = gbk_dd[bgc_id]
		# print(bgc_info)
		bgc_tid = bgc_info[0]
		if bgc_tid == 'None':
			continue
		bgc_bug = bgc_info[1]
		outaa.write('>' + 'ncbi_tid|' + str(bgc_tid) +
				'|mibig|' + bgc_id + '.1_cluster001' + '_' + m_head[1] + '_' + m_head[2] +
				'|genbank|' + m_head[6] +
				'|organism|' + bgc_bug + '\n')
		outaa.write(sequence + '\n')
	return None


def main():
	parser = make_arg_parser()
	args = parser.parse_args()
	outpath = os.path.join(args.output)
	if not os.path.isdir(outpath):
		os.mkdir(os.path.join(outpath))
		if not os.path.isdir(outpath):
			print('\nError creating output directory; check given path and try again\n')
			sys.exit()
	gbkpath = os.path.join(args.input)
	gbks = os.listdir(gbkpath)
	gbks = [f for f in gbks if f.endswith('gbk')]
	# nt = NCBITree()
	outf = open(os.path.join(outpath, 'bgc_2_tid.txt'), 'w')
	gbk_dd = defaultdict(list)
	for gbk in gbks:
		# print(gbk)
		gbk_file = os.path.join(gbkpath, gbk)
		gbk_id = gbk.split('.')[0]
		dict_list = get_tids(gbk_id, gbk_file, outf)
		gbk_dd[gbk_id] = dict_list
	outf.close()
	outaa = open(os.path.join(outpath, 'plucked_mibig_db.mpfa'), 'w')
	with open(args.mibig_aa, 'r') as inf:
		pluckify_mibig(inf, outaa, gbk_dd)
	outaa.close()


if __name__ == '__main__':
	main()
