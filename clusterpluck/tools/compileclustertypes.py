#!/usr/bin/env python
# Takes antismash output files in the current directory
# and grabs the cluster types summary file, generating
# a single txt file with all the cluster ID, types, and ranges.

# USAGE
# $ python compileclustertypes.py -o output_filename
# Returns tab-delimited txt file in the current directory.

import os, sys
import pandas as pd
import argparse


def make_arg_parser():
	parser = argparse.ArgumentParser()
	parser.add_argument("-o","--output",
		required=True,
		help="Output filename [required]") 
	return parser
def main():
	outfile = pd.DataFrame()
	for dir in os.listdir('.'):
		if dir.startswith('GCF'):
			if 'cluster_sequences' not in os.listdir(dir):
				pass
			else:
				for file in os.listdir(os.path.join(dir, 'cluster_sequences/')):
					if file.startswith('abbrev'):
						newdata = pd.read_csv(os.path.join(dir, 'cluster_sequences', file), delimiter='\t', header=0, usecols=[0,1,2])
						outfile = outfile.append(newdata)
					else:
						pass
		else:
			pass
	outfile.to_csv(args.output,sep='\t',index=False)

if __name__ == '__main__':
	parser = make_arg_parser()
	args = parser.parse_args()

	main()