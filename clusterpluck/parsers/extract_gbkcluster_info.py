#!/usr/bin/env python

# extract_gbkcluster_seq.py - version 1.0
# Extract the DNA sequence from the *.gbk output files of antiSMASH
# and return a txt file with only the concatenated sequence and FASTA-like header.

usage = 'extract_gbkcluster_info.py *.gbk'

import sys
import os
import pandas as pd
import os.path


def main():
	if len(sys.argv) < 1:  # only run if there are actually files that match
		print(usage)
		pass
	else:
		FileList= sys.argv[1:]
		Header = '>'
		#define the start of the sequence by the ORIGIN line
		seqstart = 'ORIGIN'
		sequence_begin = False
		FileNum=0
		# work through each file called by the command line
		for InfileName in FileList:
			if InfileName.endswith('final.gbk'):
				pass
			elif InfileName.endswith('.gbk'): #double check to only convert the right files
				FileNum += 1 #keep track of the number of cluster files converted
				if "cluster_sequences" not in os.listdir("."):
					os.mkdir("cluster_sequences")
				HeaderF = Header + InfileName.replace('.gbk', '')
				OutFileName = 'seq_' + InfileName + '.txt'
				OutFileName = OutFileName.replace('.gbk', '')
				OutFile = open(os.path.join('cluster_sequences', OutFileName), 'w')
				OutFile.write(HeaderF + '\n') #use the filename to ID the file on the first line
				Infile = open(InfileName, 'r')
				for line in Infile:
					if sequence_begin: #only do this if ORIGIN starts the line
						#joins together only the characters on the line in the set atcg
						OutFile.write(''.join([ch for ch in line if ch in set(('a','t','c','g'))]))
					elif line.startswith('ORIGIN'):
						sequence_begin = True #identifies the line starting with ORIGIN as sequence start
			else:
				print(usage)
				Infile.close()
				OutFile.close()
		# this loop reads the BGC summary txt document and pulls columns for ID, type, and range
		# then writes them into a new tab-delimited text file in the cluster_sequences folder
		gcfname = os.path.realpath('.')
		if "txt" in os.listdir('.'):
			for file in os.listdir('txt'):
				if file.endswith('_BGC.txt'):
					gcftmp = open('tempids.txt', 'w')
					gcfname = gcfname.split('/')[-1]
					gcfname = gcfname.replace('_genomic','')
					gcftmp.write('GCF_ID' + '\n' + gcfname)
					gcftmp.close()
					gcfdf = pd.read_csv('tempids.txt',header=0)
					BGCt = pd.read_csv(os.path.join('txt',file),delimiter='\t',header=0,usecols=[0,1,3])
					BGCtablename = "abbrev_" + file
					BGCt.insert(0,'GCF_ID',gcfdf)
					BGCof = open(os.path.join('cluster_sequences',BGCtablename),'w')
					BGCof.write(' ')
					BGCof.close()
					BGCt.to_csv(os.path.join('cluster_sequences', BGCtablename),sep='\t',index=False,na_rep=gcfname)
					os.remove('tempids.txt')
		else:
			pass
	# print to screen the number of files converted
# 	sys.stderr.write("Converted %d file(s)\n" % FileNum)

if __name__ == '__main__':
	main()
