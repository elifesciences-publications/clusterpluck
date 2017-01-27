#!/usr/bin/env Python

from ninja_utils.utils import run_command


def run_beta_div(infile, id_level, map, variable, metric, outdir, shell=False):
	"""
	Calculate beta diversity for the matched OFUs and compare between variables of interest from a mapping file.
	:param infile: the OFU table at a given cluster ID
	:param id_level: the clustering height level (0-100)
	:param map: the mapping file in standard format (#SampleID is col 1 header)
	:param variable: the mapping file header for the variable of interest
	:param metric: the distance method to use: manhattan, euclidean, canberra, bray, kulczynski, jaccard, gower, altGower, morisita, horn, mountford, raup , binomial, chao, cao, or mahalanobis.
	:param outdir: the output directory/prefix - given as a string that will be appended with the ID and .png file extension
	:param shell: whether to use the shell NOT RECOMMENDED (default=False)
	:return: the STDERR/STDOUT
	"""
	# Deal with any spaces in the file paths
	if ' ' in infile:
		infile = ''.join(['"', infile, '"'])
	if ' ' in map:
		map = ''.join(['"', map, '"'])
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
		str(map),
		str(variable),
		str(metric),
		str(outdir)]
	return run_command(cmd, shell=shell)
