#!/usr/bin/env Python

from ninja_utils.utils import run_command


def run_shogun(inpath, outpath, utree_db, threads, shell=False):
	"""
	Run SHOGUN LCA using UTree.
	:param input: Directory with the FASTA files, ".fna" extension required
	:param output: Output directory for the result files
	:param utree_db: path to the UTree database (".ctr" extension
	:param threads:
	:param cpus: the number of cpus to use
	:param shell: whether to use the shell NOT RECOMMENDED (default=False)
	:return: the STDERR/STDOUT
	"""
	# Deal with any spaces in the file paths
	if ' ' in inpath:
		inpath = ''.join(['"', inpath, '"'])
	if ' ' in outpath:
		outpath = ''.join(['"', outpath, '"'])
	if ' ' in utree_db:
		utree_db = ''.join(['"', utree_db, '"'])

	cmd = ['blastp',
		'-i', str(inpath),
		'-o', str(outpath),
		'-u', str(utree_db),
		'-p', str(threads)]
	return run_command(cmd, shell=shell)
