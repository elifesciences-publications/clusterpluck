#!/usr/bin/env Python

from ninja_utils.utils import run_command

def run_makeblastdb(infile, outfile, shell=False):
	"""
	Query your cluster sequences against an AA database (designed with MIBiG in mind) using blastp.
	:param infile: the query cluster AA sequences file (mpfa format)
	:param outfile: the resulting blastp databsefile
	:param shell: whether to use the shell NOT RECOMMENDED (default=False)
	:return: the STDERR/STDOUT
	"""
	# Deal with any spaces in the file paths
	if ' ' in infile:
		infile = ''.join(['"', infile, '"'])
	if ' ' in outfile:
		outfile = ''.join(['"', outfile, '"'])
	# DEBUG
	# print('\n')
	# print(infile)
	# print(outfile)
	# print(database)

	cmd = ['makeblastdb',
		'-in', str(infile),
		'-out', str(outfile),
		'-dbtype', 'prot',
		'-hash_index']
	# print(cmd)
	return run_command(cmd, shell=shell)
