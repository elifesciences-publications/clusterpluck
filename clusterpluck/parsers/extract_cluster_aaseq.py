#!/usr/bin/env python

# by Robin Shields-Cutler
# extract_cluster_aaseq.py - version 1.0
# August 2016
# Parses .gbk output files of antiSMASH 3.0, yielding a txt file with each module's 
# coding sequence (CDS) listed in an faa-like file. The cluster module ID is listed as:
# >NC_011593.1.cluster001_ctg1_orf00220 for example.

usage = 'extract_cluster_aaseq.py *.gbk'

import sys
import os
import os.path
import re
	
def main():
	if len(sys.argv)<1: #only run if there are actually files that match
		print(usage)
		pass
	else:
		FileList= sys.argv[1:]
		Header = '>' # start each sequence's identifier with the pacman >
		#define the start of the sequence by the CDS line
		codingstart = 'CDS   '
		title_begin = False
		sequence_begin = False
		FileNum=0
		# work through each file called by the command line
		for InfileName in FileList:
			if InfileName.endswith('final.gbk'):
				pass # avoid the 'final' summary gbk file
			elif InfileName.endswith('.gbk'): #double check to only convert the right files
				FileNum += 1 #keep track of the number of cluster files converted
				if "cluster_aa_sequences" not in os.listdir("."):
					os.mkdir("cluster_aa_sequences")
				HeaderF = Header + InfileName.replace('.gbk','')
				HeaderF = HeaderF.replace('.clu','_clu')
				OutFileName = 'aa_' + InfileName + '.mpfa'
				OutFileName = OutFileName.replace('.gbk','')
				OutFile = open(os.path.join('cluster_aa_sequences', OutFileName),'w')
				Infile = open(InfileName, 'r')
				for line in Infile:
					if title_begin: #only do this if 'CDS  ' starts the line
						if line.startswith("                     /locus_tag"):
							p = re.compile(r"^(\s+)(/locus_tag=)\"(ctg)(\d_\w+)\"")
							m = p.search(line) # searches using the regex defined above
							OutfileM = ''.join(m.group(3,4))
							OutFile.write(HeaderF + '_') #use the filename to ID the file on the first line
							OutFile.write(OutfileM + '\n')
						if line.startswith('                     /translation'):
							sequence_begin = True
						if sequence_begin:
							if line.startswith('                     /translation'):
								aa_p = re.compile(r"^(\s+)(\/translation\=\")([A-Z]+)")
								aa_m = aa_p.search(line) # searches using the regex defined above
								Outaa = aa_m.group(3)
								OutFile.write(Outaa)
							if line.startswith('                     '):
								OutFile.write(''.join([ch for ch in line if ch in set(('G,A,L,M,F,W,K,Q,E,S,P,V,I,C,Y,H,R,N,D,T'))]))
							else:
								OutFile.write('\n')
								sequence_begin = False
								if line.startswith('     CDS  '):
									title_begin = True
	# 							if line.startswith('                     /translation'):
	# 								sequence_begin = True
								else:
									title_begin = False
									sequence_begin = False
					elif line.startswith('     CDS  '):
						title_begin = True #identifies the line starting with CDS as cluster module sequence start
			else:
				print(usage)
				Infile.close()
				OutFile.close()
	
	# print to screen the number of files converted
# 	sys.stdout.write("Converted %d file(s)\n" % FileNum)
		
if __name__ == '__main__':
	main()
