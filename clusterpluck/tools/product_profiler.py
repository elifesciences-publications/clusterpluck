#!/usr/bin/env python

import argparse
import sys
import csv
import pandas as pd
import numpy as np
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
	parser.add_argument('-a', '--abundance', help='Scale the number of types by the abundance counts in the taxa table', action='store_true', required=False, default=False)
	return parser


# Make a dictionary for each organism, values are a list of the product types
def type_dict(intypes):
	type_dd = defaultdict(list)
	with open(intypes, 'r') as infile:
		reader = csv.reader(infile)
		next(reader)
		for line in reader:
			prod_type = line[1]
			cluster = line[0].split('|')[-1]
			type_dd[cluster].append(prod_type)
	# print('\nProcessed type table for %s organisms\n' % len(type_dd))
	with open(intypes, 'r') as infile2:
		types = pd.read_csv(infile2, header=0, index_col=1)
		types = set(list(types.index))
	# print('Working with %s different predicted product types\n' % len(types))
	return type_dd, types


# If matched multiple taxons (common with species annotations), just use the taxon match with the most cluster types represented
def rep_taxon(taxons, type_dd):
	longest = len(type_dd[taxons[0]])
	reptax = str(taxons[0])
	for t in taxons:
		ptypes = type_dd[t]
		if len(ptypes) > longest:
			longest = len(ptypes)
			reptax = str(t)
		else:
			continue
	return reptax


# This function returns the frequency of each product types' occurrence in each sample (not scaled by abundance counts)
def prod_types_per_sample(intaxa, type_dd, types, abundance):
	with open(intaxa, 'r') as intaxa:
		# taxa_reader = csv.reader(intaxa, delimiter='\t')
		taxa_df = pd.read_table(intaxa, delimiter='\t', header=0, index_col=0)
	n_samples = taxa_df.shape[1]
	samples = list(taxa_df.columns)
	all_types = list(types)
	compiled_df = pd.DataFrame(index=all_types)
	per_sample_taxa = []
	per_sample_strains = []
	per_sample_species = []
	per_sample_sp_st_avail = []
	for s in samples:
		sp = 0
		st = 0
		sp_st = 0
		s_name = str(s)
		ptypes = []
		# print(s_name)
		s_taxa = taxa_df[[s_name]]
		# print(s_taxa.shape[0])
		s_taxa = s_taxa.loc[(s_taxa != 0).any(axis=1)]  # Drop taxa not in this sample
		per_sample_taxa.append(s_taxa.shape[0])
		# print(s_taxa.shape[0])
		s_taxons = list(s_taxa.index)
		for taxon in s_taxons:
			taxon = str(taxon)
			if '; ' in taxon:
				';'.join(taxon.split('; '))
			# Screen through only taxa with a species or strain-level annotation
			if 's__' not in taxon or taxon.endswith('s__;t__') or taxon.endswith('s__;t__None') or taxon.endswith('s__None;t__None'):
				continue
			sp_st += 1
			# print(taxon)
			if abundance:
				taxon_abund = round(float(s_taxa.loc[taxon]))  # Save the counts for this taxon in this sample
			# If there's no strain, go to species
			if taxon.endswith('t__None') or taxon.endswith('t__'):
				k = -2
				name = taxon.split(';')[k]
				taxons = [n for n in list(type_dd.keys()) if name in n]
				if taxons:
					# print(name)
					sp += 1
			else:
				k = -1
				name = taxon.split(';')[k]
				taxons = [n for n in list(type_dd.keys()) if name in n]
				# If strain fails to match, back up and take species
				if taxons:
					# print(name)
					st += 1
				if not taxons:
					k = -2
					name = taxon.split(';')[k]
					if name.endswith('__None') or name.endswith('__') or 's__' not in name:
						continue  # If there's no species information, go to the next taxon
					taxons = [n for n in list(type_dd.keys()) if name in n]
					if taxons:
						# print(name)
						sp += 1
			# Just pick a representative taxon if multiple taxons match the species/strain name hit
			if len(taxons) > 1:
				reptax = rep_taxon(taxons, type_dd)
				prod_types = list(type_dd[reptax])
				if abundance:
					prod_types = [n for n in prod_types for _ in range(taxon_abund)]
				ptypes.extend(prod_types)
			elif len(taxons) == 1:
				prod_types = list(type_dd[str(taxons[0])])
				if abundance:  # If scaling by counts, multiply the product counts by the abundance saved earlier
					prod_types = [n for n in prod_types for _ in range(taxon_abund)]
				ptypes.extend(prod_types)
		if ptypes:
			counted = Counter(ptypes)
			samp_typedf = pd.DataFrame.from_dict(dict(counted), orient='index')
			samp_typedf.columns = [s_name]
			compiled_df = compiled_df.join(samp_typedf, how='outer')
		per_sample_species.append(sp)
		per_sample_strains.append(st)
		per_sample_sp_st_avail.append(sp_st)
	n_samples_compiled = compiled_df.shape[1]
	if n_samples != n_samples_compiled:
		print('Some samples failed to make it to the type table... started with %s, ended with %s\n' % (n_samples, n_samples_compiled))
	print('The median number of taxa per sample was %s; median number at species or strain was %s\n' % (int(np.median(per_sample_taxa)), int(np.median(per_sample_sp_st_avail))))
	print('The median number of species matches used per sample was %s; median strain matches was %s\n' % (int(np.median(per_sample_species)), int(np.median(per_sample_strains))))
	return compiled_df


def main():
	parser = make_arg_parser()
	args = parser.parse_args()
	# parse command line
	abundance = args.abundance
	intypes = args.types
	type_dd, types = type_dict(intypes)
	intaxa = args.input
	compiled_df = prod_types_per_sample(intaxa, type_dd, types, abundance)
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
