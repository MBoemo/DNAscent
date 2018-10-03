#----------------------------------------------------------
# Copyright 2017 University of Oxford
# Written by Michael A. Boemo (michael.boemo@path.ox.ac.uk)
# This software is licensed under GPL-3.0.  You should have
# received a copy of the license with this software.  If
# not, please Email the author.
#----------------------------------------------------------

import pysam
import sys
import os
import gc
import h5py
import warnings

#--------------------------------------------------------------------------------------------------------------------------------------
class arguments:
	pass


#--------------------------------------------------------------------------------------------------------------------------------------
def splashHelp():
	s = """alignHsapiens.py: Osiris preprocessing script that will align reads to the H. sapiens reference and QC on the alignment.
To run alignHsapiens.py, do:
  python alignHsapiens.py [arguments]
Example:
  python alignScerevisiae.py -r /path/to/reference.fasta --reads /path/to/reads.fastq
Required arguments are:
  -r,--reference            path to H. sapiens reference genome in fasta format,
  --reads                   path to fastq file with all reads to align.
Optional arguments are:
  -t,--threads              number of threads (default is 1 thread)."""

	print s
	exit(0)


#--------------------------------------------------------------------------------------------------------------------------------------
def parseArguments(args):

	a = arguments()
	a.threads = 1

	for i, argument in enumerate(args):
		if argument == '--reads':
			a.reads = str(args[i+1])
			
		elif argument == '-r' or argument == '--reference':
			a.reference = str(args[i+1])

		elif argument == '-t' or argument == '--threads':
			a.threads = int(args[i+1])

		elif argument == '-h' or argument == '--help':
			splashHelp()
		elif argument[0] == '-':
			splashHelp()

	#check that required arguments are met
	if not hasattr( a, 'reads') or not hasattr( a, 'reference'):
		splashHelp() 

	return a

#MAIN--------------------------------------------------------------------------------------------------------------------------------------
args = sys.argv
a = parseArguments(args)

#do the alignment with graphmap
os.system('/data/software_local/minimap2-2.10/minimap2 -ax map-ont -t ' + str(a.threads) + ' ' + a.reference + ' ' + a.reads + ' | samtools view -Sb - | samtools sort - alignments.minimap2.sorted') 
os.system('samtools index alignments.minimap2.sorted.bam')

#open the sorted bam file and output bam file
out_files = list()
sam_file = pysam.Samfile('alignments.minimap2.sorted.bam')
filtered_file = pysam.Samfile('filteredOut.minimap2.bam', "wb", template=sam_file)

#go through the sorted bam file and crop out mitochondrial and ribosomal DNA
reverseTally = 0
mDNATally = 0
rDNATally = 0
mapQTally = 0
lengthTally = 0
for record in sam_file:

	#if this read mapped
	if record.reference_id != -1:
	
		#only take reads with P(mapped to correct position) > 0.99
		if record.mapping_quality < 20:
			mapQTally += 1
			continue

		else:
			filtered_file.write( record )

sam_file.close()
filtered_file.close()

os.system('samtools index filteredOut.minimap2.bam')

sam_file = pysam.Samfile('alignments.minimap2.sorted.bam')
numOfReads = sam_file.count()
sam_file.close()

print "Total reads: ", numOfReads
print "Excluded for reverse complement: ", reverseTally, float(reverseTally)/float(numOfReads)
print "Excluded for mDNA: ", mDNATally, float(mDNATally)/float(numOfReads)
print "Excluded for rDNA: ", rDNATally, float(rDNATally)/float(numOfReads)
print "Excluded for mapping quality: ", mapQTally, float(mapQTally)/float(numOfReads)
print "Excluded for length: ", lengthTally, float(lengthTally)/float(numOfReads)
