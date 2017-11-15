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
  -d,--data                 path to top level directory of ONT reads,
  -p,--position             position of analogue in training data (valid arguments are 1and2, 3and4, or 5and6).
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

		elif argument == '-p' or argument == '--position':
			a.position = str(args[i+1])

		elif argument == '-h' or argument == '--help':
			splashHelp()
		elif argument[0] == '-':
			splashHelp()

	#check that required arguments are met
	if not hasattr( a, 'reference') or not hasattr( a, 'data') or not hasattr( a, 'position'):
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

reference = import_reference(a.reference)

if a.position == '1and2':
	analogueLoc = reference.find('NTNNNNN')
	adenineLoc = reference.find('NNNNNAN')

elif a.position == '3and4':
	analogueLoc = reference.find('NNNTNNN')
	adenineLoc = reference.find('NNNANNN')

elif a.position == '5and6':
	analogueLoc = reference.find('NNNNNTN')
	adenineLoc = reference.find('NANNNNN')

else:
	print 'Exiting with error.  Invalid argument passed to -p or --position.'
	splashHelp()

if analogueLoc == -1 or adenineLoc == -1:
	print 'Exiting with error - BrdU and/or adenine domains not found.  Did you enter the right position?'
	splashHelp()

# open an output file for each reference sequence
for x in sam_file.references:
	print x
	out_files.append(pysam.Samfile(x + ".bam", "wb", template=sam_file))

for record in sam_file:
	ref_length = sam_file.lengths[record.reference_id]

	if record.aend is None or record.query_alignment_length is None:
		continue

	ref_cover = float(record.aend - record.pos) / ref_length
	query_cover = float(record.query_alignment_length) / record.query_length

	if (record.reference_start < analogueLoc - 15) and (record.reference_end > adenineLoc + 21) and ref_cover > 0.8 and query_cover > 0.8 and (record.is_reverse == False):
		out_files[record.reference_id].write(record)


