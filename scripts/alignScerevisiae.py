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
	s = """alignScerevisiae.py: Osiris preprocessing script that will align reads to the S. cerevisiae reference and QC on the alignment.
To run alignScerevisiae.py, do:
  python alignScerevisiae.py [arguments]
Example:
  python alignScerevisiae.py -r /path/to/reference.fasta --reads /path/to/reads.fastq
Required arguments are:
  -r,--reference            path to S. cerevisiae reference genome in fasta format,
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

#do the alignment with minimap2
#os.system('/data/software_local/minimap2-2.10/minimap2 -ax map-ont -t ' + str(a.threads) + ' ' + a.reference + ' ' + a.reads + ' | samtools view -Sb - | samtools sort - alignments.minimap2.sorted') 

#do the alignment with graphmap
#os.system('graphmap align -t '+str(a.threads)+' -x sensitive -r '+a.reference+' -d ' + a.reads + ' | samtools view -Sb - | samtools sort - alignments.minimap2.sorted') 

#os.system('samtools index alignments.minimap2.sorted.bam')

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
unmappedTally = 0
numOfReads = 0
for record in sam_file:

	numOfReads += 1

	#if this read mapped
	if record.reference_id != -1:
	
		#skip mitochondrial DNA
		if record.reference_name == 'chrM':
			mDNATally += 1
			continue

		#skip ribosomal DNA
		elif record.reference_name == 'chrXII' and ( (record.reference_start > 450000 and record.reference_start < 470000) or (record.reference_end > 450000 and record.reference_end < 470000) ):
			rDNATally += 1
			continue

		#only take reads with P(mapped to correct position) > 0.99
		elif record.mapping_quality < 20:
			mapQTally += 1
			continue

		#exclude any read that's mapped to a subsequence of the reference shorter than a certain threshold
		elif record.reference_length < 20000:
			lengthTally += 1
			continue

		else:
			filtered_file.write( record )

	else:
		unmappedTally += 1

sam_file.close()
filtered_file.close()

os.system('samtools index filteredOut.minimap2.bam')

print "Total reads: ", numOfReads
print "Excluded for unmapped: ", unmappedTally, float(unmappedTally)/float(numOfReads)
print "Excluded for reverse complement: ", reverseTally, float(reverseTally)/float(numOfReads)
print "Excluded for mDNA: ", mDNATally, float(mDNATally)/float(numOfReads)
print "Excluded for rDNA: ", rDNATally, float(rDNATally)/float(numOfReads)
print "Excluded for mapping quality: ", mapQTally, float(mapQTally)/float(numOfReads)
print "Excluded for length: ", lengthTally, float(lengthTally)/float(numOfReads)

