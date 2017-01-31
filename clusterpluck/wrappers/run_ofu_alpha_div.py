#!/usr/bin/env Python

from ninja_utils.utils import run_command


def run_alpha_div(infile, id_level, mapping, variable, metric, outdir, shell=False):
	"""
	Calculate alpha diversity for the matched OFUs and compare between variables of interest from a mapping file.
	:param infile: the OFU table at a given cluster ID
	:param id_level: the clustering height level (0-100)
	:param mapping: the mapping file in standard format (#SampleID is col 1 header)
	:param variable: the mapping file header for the variable of interest
	:param metric: the diversity metric to use: shannon, simpson, or invsimpson
	:param outdir: the output directory/prefix - given as a string that will be appended with the metric, ID level, and .png file extension
	:param shell: whether to use the shell NOT RECOMMENDED (default=False)
	:return: the STDERR/STDOUT
	"""
	# Deal with any spaces in the file paths
	if ' ' in infile:
		infile = ''.join(['"', infile, '"'])
	if ' ' in mapping:
		mapping = ''.join(['"', mapping, '"'])
	if ' ' in outdir:
		outdir = ''.join(['"', outdir, '"'])
	# DEBUG
	# print('\n')
	# print(infile)
	# print(outfile)
	# print(database)

	cmd = ['ofu_beta_div.R',
		str(infile),
		str(id_level),
		str(mapping),
		str(variable),
		str(metric),
		str(outdir)]
	return run_command(cmd, shell=shell)
