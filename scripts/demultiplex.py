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
	s = """demultiplexFastq.py: Osiris preprocessing script that will demultiplex reads when basecalled to fastq.
To run demultiplexFastq.py, do:
  python demultiplexFastq.py [arguments]
Example:
  python demultiplexFastq.py -r /path/to/reference.fasta --reads /path/to/reads.fasta
Required arguments are:
  -r,--reference            path to reference file in fasta format,
  --reads                   path to fastq file with all reads to demultiplex.
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


#--------------------------------------------------------------------------------------------------------------------------------------
def displayProgress(current, total):

	barWidth = 70
	progress = float(current)/float(total)

	if progress <= 1.0:
		sys.stdout.write('[')
		pos = int(barWidth*progress)
		for i in range(barWidth):
			if i < pos:
				sys.stdout.write('=')
			elif i == pos:
				sys.stdout.write('>')
			else:
				sys.stdout.write(' ')
		sys.stdout.write('] '+str(int(progress*100))+' %\r')
		sys.stdout.flush()


#--------------------------------------------------------------------------------------------------------------------------------------
def import_reference(filename):
#	takes the filename of a fasta reference sequence and returns the reference sequence as a string.  N.B. the reference file must have only one sequence in it
#	ARGUMENTS
#       ---------
#	- filename: path to a reference fasta file
#	  type: string
#	OUTPUTS
#       -------
#	- reference: reference string
#	  type: string

	f = open(filename,'r')
	g = f.readlines()
	f.close()

	reference = ''
	for line in g:
		if line[0] != '>':
			reference += line.rstrip()
	g = None

	reference = reference.upper()

	if not all(c in ['A','T','G','C'] for c in reference):
		warnings.warn('Warning: Illegal character in reference.  Legal characters are A, T, G, and C.', Warning)

	return reference


#--------------------------------------------------------------------------------------------------------------------------------------
def split_reference(filename):

	f = open(filename,'r')
	g = f.readlines()
	f.close()	

	referenceDict = {}
	first = True

	for line in g:
		if line[0] == '>':
			if not first:
				referenceDict[key] = seq
			key = line[1:].rstrip()
			seq = ''
			first = False
		elif line == '':
			continue
		else:
			seq += line.rstrip()

	referenceDict[key] = seq		

	return referenceDict


#--------------------------------------------------------------------------------------------------------------------------------------
def print_split_reference(refDict):

	for key in refDict:
		f = open(key + '.fasta','w')
		f.write('>' + key + '\n')
		f.write(refDict[key])
		f.close()



#MAIN--------------------------------------------------------------------------------------------------------------------------------------
args = sys.argv
a = parseArguments(args)

#do the alignment with graphmap
#os.system('graphmap align -t '+str(a.threads)+' -x sensitive -r '+a.reference+' -d ' + a.reads + ' | samtools view -Sb - | samtools sort - alignments.sorted') 
#os.system('samtools index alignments.sorted.bam')

#split the reference
referenceDict = split_reference(a.reference)
print_split_reference(referenceDict)

#open an output file for each reference sequence
out_files = list()
sam_file = pysam.Samfile('alignments.sorted.bam')
for x in sam_file.references:
	print x
	out_files.append(pysam.Samfile(x + ".bam", "wb", template=sam_file))

#go through alignments.sorted.bam and sort reads into bam files corresponding to the reference they best aligned to
for record in sam_file:

	#catch bad alignment
	if record.aend is None or record.query_alignment_length is None or record.reference_name is None:
		continue

	query_cover = float(record.query_alignment_length) / float(record.query_length) #fraction of the read that aligns to the reference
	reference_cover = float( record.reference_length) / float( sam_file.lengths[record.reference_id] )
	overflow = float(record.query_length) / float( sam_file.lengths[record.reference_id] )

	#only keep reads that have the analogue ROI mapped, reads where at least 90% aligns to the reference, and reads that aren't the reverse complement
	#if query_cover > 0.9 and reference_cover > 0.9 and overflow < 1.2 and (record.is_reverse == False):

	if record.is_reverse == False:
		out_files[record.reference_id].write(record)	

#index the new BAM files
for x in sam_file.references:
	os.system('samtools index '+ x + '.bam')
