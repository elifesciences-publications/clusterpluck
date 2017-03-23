#!/usr/bin/env Python

import argparse
import pandas as pd
import csv
import os

# The arg parser
def make_arg_parser():
	parser = argparse.ArgumentParser(description='Generates an OFU table based on taxon abundance and OFU profile of those species')
	parser.add_argument('-t', '--taxons', help='The taxon table in CSV', required=True)
	parser.add_argument('-f', '--ofus', help='The OFU profile at a particular cut height, in CSV', required=True)
	parser.add_argument('-o', '--output', help='Where to save the output csv; default to "ofu_table.csv" in current working directory', required=False, default='ofu_table.csv')
	parser.add_argument('-m', '--multiples', help='When multiple strains for species-level match: '
												'[average] the OFU tallies,'
												'[summarize] all the possible OFUs,'
												'select only [universal] OFUs,'
												'select only OFUs in a [majority] of strains', required=False, default='majority')
	parser.add_argument('-p', '--profiles',
						help='Save the matching OFU profiles csv as "matching_ofu_profiles.csv" in current working directory',
						action='store_true', required=False, default=False)
	parser.add_argument('--biom', help='Also convert the final OFU table to biom format, compatible with QIIME, etc. Output must end with ".txt".', action='store_true', required=False, default=False)
	return parser


# Uses the taxon table to refine the OFU profile according to taxons that were actually found
def match_tables(taxon_file, ofu_infile, intax, ofu_index, opt):
	if taxon_file.endswith('.csv'):
		tdf = pd.read_csv(intax, header=0, index_col=0)
	if taxon_file.endswith('.txt'):
		tdf = pd.read_csv(intax, sep='\t', header=0, index_col=0)
	# odf = pd.read_csv(inofu, header=0, index_col=0)
	taxons = list(tdf.index)
	# ofu_index = list(odf.columns)
	ofu_matched = pd.DataFrame(index=ofu_index)
	# TODO: Parallelize the taxon parsing?
	for taxon in taxons:
		if len(taxon.split(';')) < 3:
			pass
		# print(taxon)
		n = []
		t_odf = pd.DataFrame(columns=ofu_index)
		taxon = str(taxon)
		if taxon.endswith('None') or taxon.endswith('__'):
			k = -2
			name = taxon.split(';')[k]
			# name = name.replace('s__', '')
		else:
			k = -1
			name = taxon.split(';')[k]
			# name = name.replace('t__', '')
		with open(ofu_infile, 'r') as inofu:
			# print(k, name)
			ofu_reader = csv.reader(inofu)
			for line in ofu_reader:
				# print(line[0])
				# print(name)
				if name in line[0]:
					# print('a match!')
					line_df = pd.DataFrame([line[1:]], columns=ofu_index, index=[line[0]], dtype='int')
					t_odf = t_odf.append(line_df)
				else:
					pass
			if t_odf.empty and len(taxon.split(';')) >= 4:
				up_name = taxon.split(';')[k - 1]
				# print(k, up_name, 'up one')
				with open(ofu_infile, 'r') as inofu:
					ofu_reader = csv.reader(inofu)
					for line in ofu_reader:
						if up_name in line[0]:
							# print('a second order match!')
							line_df = pd.DataFrame([line[1:]], columns=ofu_index, index=[line[0]], dtype='int')
							t_odf = t_odf.append(line_df)
						else:
							pass
			if t_odf.empty and len(taxon.split(';')) >= 5:
				up_name = taxon.split(';')[k - 2]
				if 'k__' or 'p__' in up_name:
					pass
				# print(k, up_name, 'up two')
				with open(ofu_infile, 'r') as inofu:
					ofu_reader = csv.reader(inofu)
					for line in ofu_reader:
						if up_name in line[0]:
							print('a third order match!')
							line_df = pd.DataFrame([line[1:]], columns=ofu_index, index=[line[0]], dtype='int')
							t_odf = t_odf.append(line_df)
						else:
							pass
			if t_odf.empty and len(taxon.split(';')) >= 6:
				up_name = taxon.split(';')[k - 3]
				if 'k__' or 'p__' in up_name:
					pass
				# print(k, up_name, 'up three')
				with open(ofu_infile, 'r') as inofu:
					ofu_reader = csv.reader(inofu)
					for line in ofu_reader:
						if up_name in line[0]:
							print('a fourth order match!')
							line_df = pd.DataFrame([line[1:]], columns=ofu_index, index=[line[0]], dtype='int')
							t_odf = t_odf.append(line_df)
						else:
							pass

		# t_odf = odf.filter(like=name, axis=0)
		if t_odf.empty:
			# print('none here')
			pass
		elif sum(t_odf.sum()) == 0:
			pass
		elif opt == 'average':  # If resolution isn't to strain, then average the OFU counts across the higher-rank group (i.e. species)
			mean = pd.DataFrame(t_odf.mean(axis=0))
			n.append(taxon)
			mean.columns = n
			ofu_matched = ofu_matched.join(mean)
		elif opt == 'summarize':  # Same as above, but here just take the summary, or the maximum number of appearances for any given OFU
			summ = pd.DataFrame(t_odf.max(axis=0))
			n.append(taxon)
			summ.columns = n
			ofu_matched = ofu_matched.join(summ)
		elif opt == 'universal':  # Same as above, but here just take the minimum set, or only OFUs shared by all strains
			univ = pd.DataFrame(t_odf.min(axis=0))
			n.append(taxon)
			univ.columns = n
			ofu_matched = ofu_matched.join(univ)
		elif opt == 'majority':  # Same as above, but here just take the majority set, or only OFUs shared by at least half of the strains
			avgm = pd.DataFrame(t_odf.mean(axis=0))
			dfbool = t_odf > 0
			maj = pd.DataFrame(dfbool.mean(axis=0))
			maj = maj >= 0.5
			maj = maj * avgm
			n.append(taxon)
			maj.columns = n
			ofu_matched = ofu_matched.join(maj)
		else:
			print('\nMust enter a valid method for dealing with multiple taxon-OFU hits, either -m "average" (default) or "summarize"\n')
			quit()
	if not ofu_matched.empty:
		ofu_matched = ofu_matched.T  # This orients the dataframe in the same way as a taxon table
		return ofu_matched
	else:
		print('\nNo OTU matches found\n')
		exit()


# Use the relative abundances to create OFU table with 'abundance' of each OFU trait
def multiply_tables(taxon_file, intaxm, ofu_matched):
	if taxon_file.endswith('.csv'):
		tdf = pd.read_csv(intaxm, header=0, index_col=0)
	if taxon_file.endswith('.txt'):
		tdf = pd.read_csv(intaxm, sep='\t', header=0, index_col=0)
	tdf = tdf.fillna(0)
	ofu_list = list(ofu_matched.index)
	if not tdf.shape[0] == ofu_matched.shape[0]:
		print('\nFYI, some taxa did not end up in the ofu profile...')
		print('taxon table dimensions =', tdf.shape[0], ',', tdf.shape[1])
		print('ofu profile dimensions = ', ofu_matched.shape[0], ',', ofu_matched.shape[1])
		print('\nUsing OFU index to limit to common taxa...\n')
		tdf = tdf.loc[ofu_list]
	print('Final taxon table dimensions =', tdf.shape[0], ',', tdf.shape[1])
	# print('Final ofu profile dimensions = ', ofu_matched.shape[0], ',', ofu_matched.shape[1], '\n')
	ofu_table = tdf.T.dot(ofu_matched)  # taking the dot product in this direction gives the relative abundance of OFUs based on observations of the taxon
	if ofu_table.empty:
		print('\nError multiplying matrices. Check dimensions.\n')
		exit()
	else:
		return ofu_table


def main():
	parser = make_arg_parser()
	args = parser.parse_args()

	# Parse command line
	with open(args.taxons, 'r') as intaxon, open(args.ofus, 'r') as inofu:
		taxon_file = args.taxons
		if taxon_file.endswith('.csv'):
			tdf = pd.read_csv(intaxon, header=0, index_col=0, usecols=[0, 1, 2])
		if taxon_file.endswith('.txt'):
			tdf = pd.read_csv(intaxon, sep='\t', header=0, index_col=0, usecols=[0, 1, 2])
		taxons = list(tdf.index)
		odf = pd.read_csv(inofu, header=0, index_col=0, nrows=2)
		ofu_index = list(odf.columns)
	opt = args.multiples
	with open(args.taxons, 'r') as intax:
		ofu_infile = args.ofus
		ofu_matched = match_tables(taxon_file, ofu_infile, intax, ofu_index, opt)
	print('Finished matching OFUs from organisms in taxon table... now multiplying for abundance...\n')
	if args.profiles:
		profile_out = 'matching_ofu_profiles.csv'
		with open(profile_out, 'w') as outf:
			ofu_matched.to_csv(outf)
	with open(args.taxons, 'r') as intaxm:
		ofu_table = multiply_tables(taxon_file, intaxm, ofu_matched)
		ofu_table = ofu_table.loc[:, (ofu_table != 0).any(axis=0)]  # removes ofus with all zeros
		print('Final ofu profile dimensions = ', ofu_table.shape[0], 'samples,', ofu_table.shape[1], ' OFUs\n')
	with open(args.output, 'w') as outf:
		ofu_table = ofu_table.round(decimals=2)
		# Make the OFU table a QIIME-compatible tsv
		ofu_table.T
		ofu_table.index.name = '#OFU ID'
		if args.output.endswith('.csv'):
			ofu_table.to_csv(outf)
		else:
			ofu_table.to_csv(outf, sep='\t')
	if args.biom and args.output.endswith('.txt'):
		print('Saving biom format file...\n')
		biom_out = '.'.join([args.output.split('.')[-2], 'biom'])
		# print(' '.join(['biom convert -i', args.output, '-o', biom_out, '--table-type="OTU table" --to-json']))
		os.system(' '.join(['biom convert -i', args.output, '-o', biom_out, '--table-type="OTU table" --to-json']))


if __name__ == '__main__':
	main()
