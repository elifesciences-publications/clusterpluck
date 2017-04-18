#!/usr/bin/env python

import argparse
import sys
import csv
import pandas as pd
from collections import defaultdict
from collections import Counter
from itertools import repeat
from multiprocessing import cpu_count
from multiprocessing import Pool


# The arg parser
def make_arg_parser():
	parser = argparse.ArgumentParser(description='Generate predicted product type table for each sample')
	parser.add_argument('-i', '--input', help='Input is pre-filtered (for low abundance taxa/low read depth) taxa table in .txt format', required=True)
	parser.add_argument('-t', '--types', help='Cluster product type table in csv', required=True)
	parser.add_argument('-o', '--output', help='Output file in tsv (.txt) or csv (.csv)', required=True)
	return parser


# Make a dictionary for each organism, values are a list of the product types
def type_dict(intypes):
	type_dd = defaultdict(list)
	with open(intypes, 'r') as infile:
		reader = csv.reader(infile)
		next(reader)
		for line in reader:
			prod_type = line[1]
			bug = line[0].split('|')[-1]
			type_dd[bug].append(prod_type)
	print('\nProcessed type table for %s organisms\n' % len(type_dd))
	with open(intypes, 'r') as infile2:
		types = pd.read_csv(infile2, header=0, index_col=1)
		types = set(list(types.index))
	print('Working with %s different predicted product types\n' % len(types))
	return type_dd, types


def prod_types_per_sample(intaxa, type_dd, types):
	with open(intaxa, 'r') as intaxa:
		# taxa_reader = csv.reader(intaxa, delimiter='\t')
		taxa_df = pd.read_table(intaxa, delimiter='\t', header=0, index_col=0)
	n_samples = taxa_df.shape[1]
	samples = list(taxa_df.columns)
	all_types = list(types)
	compiled_df = pd.DataFrame(index=all_types)
	for s in samples:
		s_name = str(s)
		ptypes = []
		# print(s_name)
		s_taxa = taxa_df[[s_name]]
		# print(s_taxa.shape[0])
		s_taxa = s_taxa.loc[(s_taxa != 0).any(axis=1)]  # Drop taxa not in this sample
		# print(s_taxa.shape[0])
		s_taxons = list(s_taxa.index)
		for taxon in s_taxons:
			taxon = str(taxon)
			if '; ' in taxon:
				';'.join(taxon.split('; '))
			# Screen through only taxa with a species or strain-level annotation
			if 's__' not in taxon and 't__' not in taxon or taxon.endswith('s__;t__') or taxon.endswith('s__;t__None'):
				continue
			# print(taxon)
			# If there's no strain, go to species
			if taxon.endswith('t__None') or taxon.endswith('t__'):
				k = -2
				name = taxon.split(';')[k]
				taxons = [n for n in list(type_dd.keys()) if name in n]
			else:
				k = -1
				name = taxon.split(';')[k]
				taxons = [n for n in list(type_dd.keys()) if name in n]
				# If strain fails to match, back up and take species
				if not taxons:
					k = -2
					name = taxon.split(';')[k]
					if name.endswith('__None') or name.endswith('__') or 's__' not in name:
						continue  # If there's no species information, go to the next taxon
					taxons = [n for n in list(type_dd.keys()) if name in n]
			for t in taxons:
				ptypes.extend(type_dd[t])
		if ptypes:
			counted = Counter(ptypes)
			samp_typedf = pd.DataFrame.from_dict(dict(counted), orient='index')
			samp_typedf.columns = [s_name]
			compiled_df = compiled_df.join(samp_typedf, how='outer')
		else:
			pass
	n_samples_compiled = compiled_df.shape[1]
	if n_samples != n_samples_compiled:
		print('Some samples failed to make it to the type table... started with %s, ended with %s\n' % (n_samples, n_samples_compiled))
	return compiled_df


def main():
	parser = make_arg_parser()
	args = parser.parse_args()
	# parse command line
	intypes = args.types
	type_dd, types = type_dict(intypes)
	intaxa = args.input
	compiled_df = prod_types_per_sample(intaxa, type_dd, types)
	compiled_df.index.name = '#TYPE ID'
	with open(args.output, 'w') as outf:
		if args.output.endswith('.txt'):
			compiled_df.to_csv(outf, sep='\t')
		elif args.output.endswith('.csv'):
			compiled_df.to_csv(outf, sep=',')
		else:
			print('Please provide either .txt or .csv output filename\n')
			sys.exit()


if __name__ == '__main__':
	main()
