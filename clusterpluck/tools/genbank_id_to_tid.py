import csv
from ninja_dojo.taxonomy import NCBITree
from clusterpluck.parsers.scrape_antismash_db import tid_to_name


def genbank_id_to_tid(gbk_id, nt_cat):
	with open(nt_cat, 'r') as nt_catalog:
		reader = csv.reader(nt_catalog, delimiter='\t')
		next(reader)
		nt = NCBITree()
		gbk_set = set()
		gbk_set.add(gbk_id)
		for line in reader:
			if line[1] in gbk_set:
				tid = line[2]
				organism = tid_to_name(tid, nt=nt)
			else:
				tid, organism = 'na', 'k__None;p__None;c__None;o__None;f__None;g__None;s__None;t__None'
	return tid, organism
