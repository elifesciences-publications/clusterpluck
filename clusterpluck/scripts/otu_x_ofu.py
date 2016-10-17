#!/usr/bin/env Python

import argparse
import sys
import pandas as pd


# The arg parser
def make_arg_parser():
	parser = argparse.ArgumentParser(description='Generates an OFU table based on taxon abundance and OFU profile of those species')
	parser.add_argument('-t', '--taxons', help='The taxon table in CSV', required=True)
	parser.add_argument('-f', '--ofus', help='The OFU profile at a particular cut height, in CSV', required=True)
	parser.add_argument('-o', '--output', help='Where to save the output csv; default to "ofu_table.csv" in current working directory', required=False, default='ofu_table.csv')
	parser.add_argument('-m', '--multiples', help='When multiple strains for species: [average] the OFU tallies or [summarize] all the possible OFUs', required=False, default='average')
	return parser


# Uses the taxon table to refine the OFU profile according to taxons that were actually found by SHOGUN
def match_tables(intax, inofu, opt):
	tdf = pd.read_csv(intax, header=0, index_col=0)
	odf = pd.read_csv(inofu, header=0, index_col=0)
	taxons = list(tdf.index)
	ofu_index = list(odf.columns)
	ofu_matched = pd.DataFrame(index=ofu_index)
	for taxon in taxons:
		n = []
		taxon = str(taxon)
		if taxon.endswith('t__'):
			name = taxon.split(';')[-2]
			name = name.replace('s__', '')
		else:
			name = taxon.split(';')[-1]
			name = name.replace('t__', '')
		t_odf = odf.filter(like=name, axis=0)
		if t_odf.empty:
			pass
		elif opt == 'average':  # If resolution isn't to strain, then average the OFU counts across the higher-rank group (i.e. species)
			mean = pd.DataFrame(t_odf.mean(axis=0))
			n.append(taxon)
			mean.columns = n
			ofu_matched = ofu_matched.join(mean)
		elif opt == 'summarize':  # Same as above, but here instead of averaging, just take the summary, or the maximum number of appearances for any given OFU
			summ = pd.DataFrame(t_odf.max(axis=0))
			n.append(taxon)
			summ.columns = n
			ofu_matched = ofu_matched.join(summ)
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
def multiply_tables(intaxm, ofu_matched):
	tdf = pd.read_csv(intaxm, header=0, index_col=0)
	tdf = tdf.fillna(0)
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
	with open(args.taxons, 'r') as intax:
		with open(args.ofus, 'r') as inofu:
			opt = args.multiples
			ofu_matched = match_tables(intax, inofu, opt)
		with open(args.taxons, 'r') as intaxm:
			ofu_table = multiply_tables(intaxm, ofu_matched)
		with open('matching_ofu_profiles.csv', 'w') as profile_out:
			ofu_matched.to_csv(profile_out)
		with open(args.output, 'w') as outf:
			ofu_table.to_csv(outf)


if __name__ == '__main__':
	main()
