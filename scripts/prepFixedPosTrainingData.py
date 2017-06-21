#----------------------------------------------------------
# Copyright 2017 University of Oxford
# Written by Michael A. Boemo (michael.boemo@path.ox.ac.uk)
#----------------------------------------------------------


import numpy as np
import sys
import warnings
import h5py
import pysam
import re
import os
import gc
import math
from joblib import Parallel, delayed #for parallel processing
import multiprocessing #for parallel processing


#--------------------------------------------------------------------------------------------------------------------------------------
class arguments:
	pass


#--------------------------------------------------------------------------------------------------------------------------------------
def splashHelp():
	s = """prepFixedPosTrainingData.py: Osiris preprocessing script that will format ONT reads in the Osiris training format.
To run prepFixedPosTrainingData.py, do:
  python prepFixedPosTrainingData.py [arguments]
Example:
  python prepFixedPosTrainingData.py -r /path/to/reference.fasta -p 45 -m /path/to/template_median68pA.5mer.model -d /path/to/reads -o output.foh -t 20
Required arguments are:
  -r,--reference            path to reference file in fasta format,
  -p,--position             position of analogue in the reference (indexing starts from 0),
  -d,--data                 path to top level directory of ONT reads,
  -m,--5mer-model           path to 5mer pore model file (provided by ONT) to normalise reads,
  -o,--output               path to the output training .foh or detection .fdh file.
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
			a.position = int(args[i+1])

		elif argument == '-m' or argument == '--5mer-model':
			a.fiveMerModel = str(args[i+1])

		elif argument == '-o' or argument == '--output':
			a.outFoh = str(args[i+1])

		elif argument == '-h' or argument == '--help':
			splashHelp()

	return a


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

	#empty reads string, and count the number of subdirectories so we can print progress
	reads = ''
	numSubdirectories = len(next(os.walk(pathToReads, topdown=True))[1])
	readCount = 0

	#recursively go through the directory and subdirectories and extract fasta seqs until you reach the buffer, then write, release, and garbage collect
	for root, dirs, files in os.walk(pathToReads, topdown=True):

		for fast5file in files:

			readCount += 1

			if fast5file.endswith('.fast5'):
		
				#print progress every 5 subdirectories of reads
				if readCount % 10000 == 0:
					sys.stdout.write("\rExporting fast5 reads to fasta... read " + str(readCount))
					sys.stdout.flush()

				try:
					#open the fast5 file with h5py and grab the fastq
					ffast5 = h5py.File(root+'/'+fast5file,'r')
					fastq = ffast5[fast5path2fastq].value
					ffast5.close()
					fasta = fastq.split('\n')[1]
			
					#append the sequence in the fasta format, with the full path to the fast5 file as the sequence name
					reads += '>'+root+'/'+fast5file+'\n'+fasta+'\n'

				except KeyError:
					warnings.warn('File '+root+'/'+fast5file+' did not have a valid fastq path.  Skipping.', Warning)

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


#--------------------------------------------------------------------------------------------------------------------------------------
def export_trainingDataToFoh( kmer2normalisedReads, filename ):
#	takes reads from a run, aligns them to a reference, and separates the resulting bam file by each reference
#	ARGUMENTS
#       ---------
#	- kmer2normalisedReads: output of import_FixedPosTrainingData
#	  type: list of lists
#	OUTPUTS
#       -------
#	- exports read observation data in the .foh format, to be read by C++ Osiris
	
	f = open(filename,'w')
	f.write( '>FixedPositionTrainingData\n' )

	for read in kmer2normalisedReads[key]:
			
		readsStr = map( str, read )

		f.write( ' '.join(readsStr) + '\n' )

	f.close()


#--------------------------------------------------------------------------------------------------------------------------------------
def serial_calculate_normalisedEvents(fast5Files, poreModelFile):
#	For fixed position training data - small enough to be done in serial.  Uses a 5mer model to calculate the shift and scale
#	pore-specific parameters for each individual read in a list of fast5 files
#	ARGUMENTS
#       ---------
#	- fast5Files: list of fast5 files whose events should be normalised
#	  type: list of strings
#	- poreModelFile: path to a pore model file.  This should be a 5mer model in the ONT format
#	  type: string
#	OUTPUTS
#       -------
#	- allNormalisedReads: a list, where each member is itself a list of events that have been normalised to the pore model
#	  type: list

	#open the 5mer model and make a map from the 5mer (string) to [mean,std] (list)
	kmer2MeanStd = import_poreModel(poreModelFile)

	#now iterate through all the relevant fast5 files so we only need to open the model file once
	allNormalisedReads = []
	for f5File in fast5Files:

		#use 5mer model to calculate shift, scale, drift, and var to normalise events for the pore
		f = h5py.File(f5File,'r')
		path = '/Analyses/Basecall_1D_000/BaseCalled_template/Events'
		Events = f[path]
		A = np.zeros((2,2))
		b = np.zeros((2,1))
		for event in Events:
			 if float(event[7]) > 0.30: #if there's a high probability (>30%) that the 5mer model called by Metrichor was the correct one
				model_5mer = event[4]
				event_mean = float(event[0])
				model_mean = kmer2MeanStd[model_5mer][0]
				model_std = kmer2MeanStd[model_5mer][1]
				
				#update matrix A
				A[0,0] += 1/(model_std**2)
				A[1,0] += 1/(model_std**2)*model_mean
				A[1,1] += 1/(model_std**2)*model_mean**2

				#update vector b
				b[0] += 1/(model_std**2)*event_mean
				b[1] += 1/(model_std**2)*event_mean*model_mean

		#use symmetry of A
		A[0,1] = A[1,0]

		#solve Ax = b to find shift and scale
		x = np.linalg.solve(A,b)
		shift = x[0][0]
		scale = x[1][0]

		#go through the same events as before and normalise them to the pore model using scale and shift
		normalisedEvents = []
		for event in Events:
			if float(event[7]) > 0.30: #if there's a high probability (>30%) that the 5mer model called by Metrichor was the correct one
				event_mean = float(event[0])
				normalisedEvents.append( event_mean/scale - shift)

		allNormalisedReads.append(normalisedEvents)

		f.close()
	
	return allNormalisedReads


#--------------------------------------------------------------------------------------------------------------------------------------
def import_FixedPosTrainingData(bamFile, poreModelFile):
#	Used to import training data from reads that have an analogue in a fixed context.
#	Creates a map from kmer (string) to a list of lists, where each list is comprised of events from a read
#	First reads a BAM file to see which reads (readIDs, sequences) aligned to the references based on barcoding.  Then finds the fast5 files
#	that they came from, normalises the events to a pore model, and returns the list of normalised events.
#	ARGUMENTS
#       ---------
#	- bamFile: a BAM file from the alignment
#	  type: string
#	- poreModelFile: ONT model file for 5mers that can be used to normalised for shift and scale
#	  type: string
#	OUTPUTS
#       -------
#	- normalisedReads: a list of lists, where each element is a list of normalised events for a given read
#	  type: list

	#open up the BAM file that has been sorted by the reference that we're interested in
	f = pysam.AlignmentFile(bamFile,'r')

	#count the records in the bam file
	numRecords = f.count()
	print str(numRecords) + ' records in BAM file.'

	#iterate through the bam file, and for every record, add the path to the fast5 file to a list
	fast5files = []
	f = pysam.AlignmentFile(bamFile,'r')
	for record in f:

		fast5files.append(record.query_name)

	f.close()

	#hand this list of fast5 files to calculate_normalisedEvents which will normalise them for shift and scale
	normalisedReads = serial_calculate_normalisedEvents(fast5files, poreModelFile)

	return normalisedReads


#--------------------------------------------------------------------------------------------------------------------------------------
def alignAndSort(readsDirectory, pathToReference, position, threads=1):
#	takes reads from a run, aligns them to a reference, and separates the resulting bam file by each reference
#	ARGUMENTS
#       ---------
#	- readsDirectory: full path to the directory that contains the fast5 files for the run 
#	  type: string
#	- pathToReference: full path to the reference file that has the reference (or references, if they're barcoded) for all sequences present in the run
#	  type: string
#	- position: where the analogue is in the 7mer (1&2, 3&4, 5&6)
#	  type: string
#	- threads: number of threads on which to run BWA-MEM 
#	  type: int 
#	OUTPUTS
#       -------
#	- (BAMrecords, redundant_A_Loc): tuple consisting of a list of BAM records and the location of the redundant adenine
#	  type: tuple of (list, int)
	
	#take the directory with the fast5 reads in it and export them to a fasta file in the current working directory
	import_fasta(readsDirectory, os.getcwd()+'/reads.fasta')

	#index the reference
	os.system('bwa index ' + pathToReference)

	#align the reads.fasta file created above to the reference with bwa-mem, then sort the bam file
	os.system('bwa mem -t '+str(threads)+' -x ont2d '+pathToReference+' reads.fasta | samtools view -Sb - | samtools sort - alignments.sorted') 
	os.system('samtools index alignments.sorted.bam')

	sam_file = pysam.Samfile('alignments.sorted.bam')
	reference = import_reference(pathToReference)
	BAMrecords = []

	for record in sam_file:

		ref_length = sam_file.lengths[record.reference_id]

		if (record.aend is None) or (record.query_alignment_length is None) or record.is_reverse:
			continue

		#if the alignment indicates that this read doesn't have the critical regions, then ignore it and don't use it for training		
		if (record.reference_start > position - 10) or (record.reference_end < position + 10) or (record.query_alignment_length/ref_length < 0.8):
			continue
			
		BAMrecords.append(record)

	return BAMrecords


#MAIN--------------------------------------------------------------------------------------------------------------------------------------
#parse arguments
args = sys.argv
a = parseArguments(args)

#do alignment QC
BAMrecords = alignAndSort(a.data, a.reference, a.position, a.threads)

#normalise the training data according to the ONT 5mer model
trainingData = import_FixedPosTrainingData(BAMrecords, a.fiveMerModel)

#write normalised reads to file in the .foh format
export_trainingDataToFoh( trainingData, a.outFoh )

