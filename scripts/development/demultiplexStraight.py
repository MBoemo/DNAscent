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
	s = """demultiplexHairpin.py: Osiris preprocessing script that will demultiplex barcoded reads for hairpin training.
To run demultiplexHairpin.py, do:
  python demultiplexHairpin.py [arguments]
Example:
  python demultiplexHairpin.py -r /path/to/reference.fasta -d /path/to/reads -p 34 -t 20
Required arguments are:
  -r,--reference            path to reference file in fasta format,
  -d,--data                 path to top level directory of ONT reads.
Optional arguments are:
  -t,--threads              number of threads (default is 1 thread)."""

	print s
	exit(0)


#--------------------------------------------------------------------------------------------------------------------------------------
def parseArguments(args):

	a = arguments()

	for i, argument in enumerate(args):
		if argument == '-t' or argument == '--threads':
			a.threads = int(args[i+1])
			
		elif argument == '-r' or argument == '--reference':
			a.reference = str(args[i+1])

		elif argument == '-d' or argument == '--data':
			a.data = str(args[i+1])

		elif argument == '-h' or argument == '--help':
			splashHelp()
		elif argument[0] == '-':
			splashHelp()

	#check that required arguments are met
	if not hasattr( a, 'reference') or not hasattr( a, 'data'):
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

	if not all(c in ['A','T','G','C','N'] for c in reference):
		warnings.warn('Warning: Illegal character in reference.  Legal characters are A, T, G, C, and N.', Warning)

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


#--------------------------------------------------------------------------------------------------------------------------------------
def import_fasta(pathToReads, outFastaFilename):
#	takes a directory with fast5 nanopore reads at the top level, and extracts the 2D sequences in fasta format with the path to the file as the fasta header
#	ARGUMENTS
#       ---------
#	- pathToReads: full path to the directory that contains the fast5 files
#	  type: string
#	- outFastaFilename: filename for the output fasta file that contains all of the reads
#	  type: string
#	OUTPUTS
#       -------
#	- a fasta file written to the directory specified

	buffersize = 1024

	#output file to write on
	fout = open(outFastaFilename,'w')

	#path through the fast5 tree to get to the fastq sequence
	fast5path2fastq = '/Analyses/Basecall_1D_000/BaseCalled_template/Fastq'
	fast5path2score = '/Analyses/Basecall_1D_000/Summary/basecall_1d_template'

	#empty reads string, and count the number of subdirectories so we can print progress
	reads = ''
	numSubdirectories = len(next(os.walk(pathToReads, topdown=True))[1])
	readCount = 0
	failedCounter = 0

	#recursively go through the directory and subdirectories and extract fasta seqs until you reach the buffer, then write, release, and garbage collect
	for root, dirs, files in os.walk(pathToReads, topdown=True):

		for fast5file in files:

			readCount += 1

			if fast5file.endswith('.fast5'):
		
				#print progress every 5 subdirectories of reads
				if readCount % 10000 == 0:
					sys.stdout.write("\rExporting fast5 reads to fasta... read " + str(readCount) + ' with '+str(failedCounter)+' failed to open.')
					sys.stdout.flush()

				try:
					#open the fast5 file with h5py and grab the fastq
					ffast5 = h5py.File(root+'/'+fast5file,'r')
					fastq = ffast5[fast5path2fastq].value
					score = float(ffast5[fast5path2score].attrs.__getitem__('mean_qscore'))
					ffast5.close()
					fasta = fastq.split('\n')[1]
			
					#append the sequence in the fasta format, with the full path to the fast5 file as the sequence name
					#uncomment below to threshold by read quality score
					#if score > 10:
					reads += '>'+root+'/'+fast5file+'\n'+fasta+'\n'

				except KeyError:
					#warnings.warn('File '+root+'/'+fast5file+' did not have a valid fastq path.  Skipping.', Warning)
					failedCounter += 1
					pass

				except IOError:
					warnings.warn('File '+root+'/'+fast5file+' could not be opened and may be corrupted.  Skipping.', Warning)

				#write to the file and release the buffer
				if readCount % buffersize == 0:
					fout.write(reads)
					fout.flush()
					os.fsync(fout .fileno())
					reads = ''
					gc.collect()

		#flush the buffer and write once we're reached the end of fast5 files in the subdirectory
		fout.write(reads)
		fout.flush()
		os.fsync(fout .fileno())
		reads = ''
		gc.collect()
	
	#close output fasta file	
	fout.close()
	print 'Could not open ' + str(failedCounter) + ' reads.'


#MAIN--------------------------------------------------------------------------------------------------------------------------------------
args = sys.argv
a = parseArguments(args)

import_fasta(a.data, os.getcwd()+'/reads.fasta')

os.system('bwa index ' + a.reference)
os.system('graphmap align -t '+str(a.threads)+' -x sensitive -r '+a.reference+' -d reads.fasta | samtools view -Sb - | samtools sort - alignments.sorted') 
os.system('samtools index alignments.sorted.bam')


sam_file = pysam.Samfile('alignments.sorted.bam')
out_files = list()

referenceDict = split_reference(a.reference)
print_split_reference(referenceDict)
posDict = {}

for key in referenceDict:

	if referenceDict[key].find('NTNNNNN') != -1:
		posDict[key] = referenceDict[key].find('NTNNNNN')
	elif referenceDict[key].find('NNNTNNN') != -1:
		posDict[key] = referenceDict[key].find('NNNTNNN')
	elif referenceDict[key].find('NNNNNTN') != -1:
		posDict[key] = referenceDict[key].find('NNNNNTN')
	else:
		print 'Exiting with error - BrdU and/or adenine domains not found.  Is the right domain in the reference?'
		splashHelp()

# open an output file for each reference sequence
for x in sam_file.references:
	print x
	out_files.append(pysam.Samfile(x + ".bam", "wb", template=sam_file))

for record in sam_file:
	ref_length = sam_file.lengths[record.reference_id]

	if record.aend is None or record.query_alignment_length is None or record.reference_name is None:
		continue

	analogueLoc = posDict[record.reference_name]

	query_cover = float(record.query_alignment_length) / record.query_length

	#only keep reads that have the analogue ROI mapped, reads where at least 80% aligns to the reference, and reads that aren't the reverse complement
	if (record.reference_start < analogueLoc - 15) and (record.reference_end > analogueLoc + 21) and query_cover > 0.8 and (record.is_reverse == False):
		out_files[record.reference_id].write(record)

#index the new BAM files
for x in sam_file.references:
	os.system('samtools index '+ x + '.bam')
