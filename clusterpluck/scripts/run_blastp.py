#!/usr/bin/env Python

from ninja_utils.utils import run_command


def run_blastp(infile, outfile, database, cpus, evalue=1e-10, shell=False):
	"""
	Query your cluster sequences against an AA database (designed with MIBiG in mind) using blastp.
	:param infile: the query cluster AA sequences file (mpfa format)
	:param outfile: the resulting blastp output table file
	:param database: path to the blast database files, generated with makeblastdb command
	:param evalue: the evalue cutoff (default=1e-10)
	:param cpus: the number of cpus to use
	:param shell: whether to use the shell NOT RECOMMENDED (default=False)
	:return: the STDERR/STDOUT
	"""
	cmd = ['blastp',
		'-db', str(database),
		'-query', str(infile),
		'-out', str(outfile),
		'-num_threads', str(cpus),
		'-outfmt', '6',
		'-max_hsps', '1',
		'-evalue', evalue]
	return run_command(cmd, shell=shell)
