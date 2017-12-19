#!/usr/bin/env Python

import argparse
import os
import csv
import datetime
from shutil import rmtree
import numpy as np
from collections import defaultdict
from multiprocessing import cpu_count


# Arg parser
def make_arg_parser():
	parser = argparse.ArgumentParser(description='Given the processed shotgun reads, measure coverage of predicted BGCs in metagenome')
	parser.add_argument('-s', '--seqs', help="The combined seqs file (fasta) containing QC'd shotgun reads", required=False)
	parser.add_argument('--bgc_fna', help="The BGC DNA fasta file (all BGCs in one multi-fasta file)", required=False)
	parser.add_argument('--bgc_edx', help="The pre-made BURST database file)", required=False)
	parser.add_argument('--bgc_acx', help="The pre-made BURST accelerator file)", required=False)
	parser.add_argument('-a', '--alignment', help='A pre-calculated all-vs-all alignment file in .b6 format', required=False)
	parser.add_argument('--debug', help='Keep the temp directory for debugging', default=False, action='store_true', required=False)
	# parser.add_argument('-b', '--bread', help='Where to find the header for the sequence (default="ref|,|")', default='ref|,|')
	parser.add_argument('-l', '--lengths', help="The table providing length of each BGC in bp", required=True)
	parser.add_argument('-o', '--output', help='Directory in which to save the results (default = cwd)', required=False, default='.')
	# parser.add_argument('-t', '--threshold', help='Coverage threshold at which to accept a cluster alignment (default = 75%)', required=False, default=75, type=float)
	parser.add_argument('-p', '--threads', help='Number of processors to use for alignment (default = all available)', required=False, default='-')
	return parser


# Make a database if necessary
def create_bgc_db(bgc_db, temp_path, burst):
	acc = os.path.join(temp_path, 'bgc_db.acc')
	edb = os.path.join(temp_path, 'bgc_db.edb')
	os.system(' '.join([burst, '-r', bgc_db, '-a', acc, '-o', edb, '-d DNA -s']))
	os.system(' '.join([burst, '-a', acc, '-o /dev/null']))
	os.system(' '.join([burst, '-r', edb, '-o /dev/null']))
	os.remove(acc)
	os.remove(edb)
	acx = os.path.join(temp_path, 'bgc_db.acx')
	edx = os.path.join(temp_path, 'bgc_db.edx')
	if not os.path.isfile(acx) or not os.path.isfile(edx):
		raise ValueError('BURST databases not successfully generated: check input fasta and BURST installation')
	return acx, edx


# Run the BURST alignment
def align_bgcs(temp_path, seqs, acx, edx, cpus):
	alignment = os.path.join(temp_path, 'bgc_alignments.b6')
	if acx != '-':
		os.system(' '.join(['burst15 -m ALLPATHS -fr -hr -q', seqs, '-n -r', edx, '-o', alignment, '-t', cpus, '-i 0.9', '-a', acx]))
	else:
		os.system(' '.join(['burst15 -m ALLPATHS -fr -hr -q', seqs, '-n -r', edx, '-o', alignment, '-t', cpus, '-i 0.9']))
	print("BURST alignment completed")
	return alignment


def get_percent_coverage(tally):
	# Get number of hits
	hits = np.sum(tally > 0)
	# Get number of misses
	misses = np.sum(tally == 0)
	# Coverage
	percent_coverage = hits/(hits+misses)*100
	return percent_coverage


def zero_runs(tally):
	# Stack Overflow:
	# https://stackoverflow.com/questions/24885092/finding-the-consecutive-zeros-in-a-numpy-array
	# Create an array that is 1 where a is 0, and pad each end with an extra 0.
	iszero = np.concatenate(([0], np.equal(tally, 0).view(np.int8), [0]))
	absdiff = np.abs(np.diff(iszero))
	# Runs start and end where absdiff is 1.
	ranges = np.where(absdiff == 1)[0].reshape(-1, 2)
	return ranges


def max_uncovered_region(tally):
	coverages = zero_runs(tally)
	if coverages.size == 0:
		return 0
	else:
		return np.max(coverages[:, 1] - coverages[:, 0])


def expected_coverage(tally, median_read_length):
	number_of_reads = np.sum(tally)/median_read_length
	return 1-np.exp(-(number_of_reads*median_read_length)/tally.shape[0])


def invert_bincount(c):
	return np.repeat(np.arange(c.size), c)

# plt.hist(invert_bincount(sample_dict[gene_cluster_id]), bins=100)


# Run across all samples to get coverage in the meta-metagenome
def meta_coverage(sample_to_gene_cluster_to_row):
	gene_cluster_to_row = dict()
	for key, value in sample_to_gene_cluster_to_row.items():
		for gene_cluster, row in sample_to_gene_cluster_to_row[key].items():
			if gene_cluster not in gene_cluster_to_row:
				gene_cluster_to_row[gene_cluster] = row
			else:
				gene_cluster_to_row[gene_cluster] += row

# tally = gene_cluster_to_row['ncbi_tid|451515|genbank|CP000255.1_cluster003|organism|k__Bacteria;p__Firmicutes;c__Bacilli;o__Bacillales;f__Staphylococcaceae;g__Staphylococcus;s__Staphylococcus_aureus;t__Staphylococcus_aureus_subsp._aureus_USA300_FPR3757']


def main():
	parser = make_arg_parser()
	args = parser.parse_args()
	outpath = args.output
	if args.seqs:
		seqs = args.seqs
	lengths = args.lengths
	if args.threads == '-':
		cpus = int(cpu_count())
	else:
		cpus = int(args.threads)
	# if not args.metrics == '-':
	# 	output_metrics = str(args.metrics)
	# else:
	# 	args.metrics = 'all'
	# threshold = float(args.threshold)
	if not os.path.isdir(outpath):
		os.mkdir(outpath)
		if not os.path.isdir(outpath):
			raise ValueError('Error creating output directory; check given path and try again')
	timestamp = datetime.datetime.now()
	timestamp = timestamp.strftime("%Y%m%d_%H%M")
	tempdir = '_'.join(['tempz', timestamp])
	temp_path = os.path.join(outpath, tempdir)
	os.mkdir(temp_path)
	if args.bgc_fna:
		fna_bytes = int(os.path.getsize(args.bgc_fna))
		if fna_bytes < 1000000000:
			burst = 'burst12'
		else:
			burst = 'burst15'
		bgc_db = args.bgc_fna
		acx, edx = create_bgc_db(bgc_db, temp_path, burst)
	elif args.bgc_edx:
		edx = args.bgc_edx
		if args.bgc_acx:
			acx = args.bgc_acx
		else:
			acx = '-'
	elif args.alignment:
		alignment = args.alignment
	elif not args.bgc_fna and not args.bgc_edx and not args.alignment:
		raise ValueError('No reference sequences or database provided.\nPlease include something to align to or the pre-computed alignment file')
	# Generate the alignment between sequences and BGCs
	if not alignment:
		try:
			os.path.isfile(seqs)
		except ValueError:
			print('Sequences file not found')
		alignment = align_bgcs(temp_path, seqs, acx, edx, cpus)

	# Calculate coverage based on the alignment #
	# Parse the BGC lengths file
	gene_cluster_sizes = defaultdict(int)
	with open(lengths) as inf:
		csv_inf = csv.reader(inf, delimiter='\t')
		for row in csv_inf:
			gene_cluster_sizes[row[0]] = int(row[-1])

	# Parse the alignment file
	sample_to_gene_cluster_to_row = defaultdict(dict)
	sample_names = []
	read_lengths = []
	sample_cluster_hit_counter = defaultdict(dict)
	with open(alignment) as inf:
		csv_inf = csv.reader(inf, delimiter="\t")
		for line in csv_inf:
			sample_id = line[0].split('_')[0]
			sample_names.append(sample_id)
			sample_dict = sample_to_gene_cluster_to_row[sample_id]
			bgc_hit_counter = sample_cluster_hit_counter[sample_id]
			gene_cluster_id = line[1]
			if not gene_cluster_id in sample_dict:
				sample_dict[gene_cluster_id] = np.zeros(gene_cluster_sizes[gene_cluster_id], dtype=int)
				bgc_hit_counter[gene_cluster_id] = 0
			begin = int(line[9])
			end = int(line[8])
			# Reverse complement
			if begin > end:
				begin = int(line[8])
				end = int(line[9])
			read_lengths.append(end-begin)
			sample_dict[gene_cluster_id][begin:end] += 1
			bgc_hit_counter[gene_cluster_id] += 1
	read_lengths_array = np.array(read_lengths, dtype=int)
	median_read_length = np.median(read_lengths_array)

	# Loop through the data to generate the output table
	# coverage_result = defaultdict(dict)
	outfile = os.path.join(outpath, 'coverage_result_%s.txt') % timestamp
	with open(outfile, 'w') as outf:
		outf.write('sample_id\tbgc_id\tbgc_length\tmax_gap\thits_to_bgc\ttotal_base_coverage\tpercent_coverage\texpected_coverage\tratio_covered_to_expected\n')
		sample_names = set(sample_names)
		for sample in sample_names:
			sample_dict_parse = sample_to_gene_cluster_to_row[sample]
			sample_hit_counts = sample_cluster_hit_counter[sample]
			for bgc_id in sample_dict_parse:
				tally = sample_dict_parse[bgc_id]
				hits_to_bgc = sample_hit_counts[bgc_id]
				total_base_coverage = int(np.sum(tally))  # Need this?
				max_gap = max_uncovered_region(tally)
				# print(max_gap)
				percent_coverage = get_percent_coverage(tally)
				# print(percent_coverage)
				predicted_coverage = expected_coverage(tally, median_read_length) * 100
				# print(predicted_coverage)
				coverage_ratio = percent_coverage / predicted_coverage
				# print(coverage_ratio)
				outf.write('%s\t%s\t%d\t%d\t%d\t%d\t%s\t%s\t%s\n' % (sample, bgc_id, gene_cluster_sizes[bgc_id], max_gap, hits_to_bgc, total_base_coverage, percent_coverage, predicted_coverage, coverage_ratio))
	print('Coverage analysis complete. Find your results at:\n%s' % outfile)
	if not args.debug:
		rmtree(temp_path)
	# The calls to get the values #
	# max_gap = max_uncovered_region(sample_dict[gene_cluster_id])
	# percent_coverage = get_percent_coverage(sample_dict[gene_cluster_id])
	# expected_coverage = expected_coverage(sample_dict[gene_cluster_id], median_read_length) * 100
	# coverage_ratio = percent_coverage / expected_coverage


if __name__ == '__main__':
	main()
