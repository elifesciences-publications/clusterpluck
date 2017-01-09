#!/usr/bin/env Python

import argparse
import sys
import pandas as pd
import scipy.spatial.distance as ssd
from collections import defaultdict
from scipy.cluster.hierarchy import linkage
from scipy.cluster.hierarchy import cut_tree


# The arg parser
def make_arg_parser():
	parser = argparse.ArgumentParser(description='Performs hierarchical clustering on the percent identity BGC scores matrix')
	parser.add_argument('-i', '--input', help='Input file: A percent identity (pident) scores matrix of strain_cluster comparisons', required=True)
	parser.add_argument('-o', '--output', help='Where to save the output csv; default to screen', required=False, default='-')
	parser.add_argument('-t', '--height', help='At what height to cut the tree, default is 0.3', required=False, default=0.3)
	parser.add_argument('-c', '--ofu', help='Comma-separated list of the ofus (e.g. ofu00001,ofu00003) to reveal', required=False, type=str)
	return parser


# Perform the hierarchical clustering
def process_hierarchy(inf, h):
	df = pd.read_csv(inf, header=0, index_col=0)
	df = df.fillna(0)
	strains = df.index
	df = 1 - (df / 100)
	df_v = ssd.squareform(df, force='tovector', checks=False)  # flatten matrix to condensed distance vector
	li = linkage(df_v)  # create the dendrogram from the condensed distance vector
	hclus = cut_tree(li, height=h)  # using the height (percent ID as decimal, for example), cluster OFUs from dendrogram
	hclus = pd.DataFrame(hclus, index=strains)
	hclus.ix[:, 0] += 1  # cut_tree defaults to the first 'cluster' being named "0"; this just bumps all IDs +1
	return hclus


def bgcs_in_ofu(ofus, hclus):
	dd = defaultdict(list)
	for value, key in hclus.itertuples(index=True):
		key = str('%05d' % key)
		dd[key].extend(list([value]))
	ofu_list = ofus.split(',')
	for ofu in ofu_list:
		ofu = str(ofu)
		if ofu.startswith('ofu'):
			ofu_n = str(ofu.replace('ofu', ''))
		bgcs = dd[ofu_n]
		print('\nThis ofu cluster, %s, contains the following BGCs:' % ofu)
		print(bgcs, '\n')
	return None


def main():
	parser = make_arg_parser()
	args = parser.parse_args()
	# Parse command line
	with open(args.input, 'r') as inf:
		h = float(args.height)
		hclus = process_hierarchy(inf, h)
	if args.ofu:
		ofus = args.ofu
		bgcs_in_ofu(ofus, hclus)
		sys.exit()
	with open(args.output, 'w') if args.output != '-' else sys.stdout as outf:
		hclus.to_csv(outf)

if __name__ == '__main__':
	main()
