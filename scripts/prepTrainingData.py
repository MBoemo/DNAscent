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
	s = """prepTrainingData.py: Osiris preprocessing script that will format ONT reads in the Osiris training/detection format.
To run prepTrainingData.py, do:
  python prepData.py [arguments]
Example:
  python prepTrainingData.py -r /path/to/reference.fasta -p 3&4 -m /path/to/template_median68pA.5mer.model -d /path/to/reads -o output.foh -t 20
Required arguments are:
  -r,--reference            path to reference file in fasta format,
  -p,--position             position of analogue in training data (valid arguments are 12, 34, or 56),
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
			a.position = str(args[i+1])

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
#	- kmer2normalisedReads: output of import_HairpinTrainingData 
#	  type: dictionary
#	OUTPUTS
#       -------
#	- exports read observation data in the .foh format, to be read by C++ Osiris
	
	f = open(filename,'w')

	for key in kmer2normalisedReads:
		
		f.write( '>'+key+'\n' )

		for read in kmer2normalisedReads[key]:
			
			readsStr = map( str, read )

			f.write( ' '.join(readsStr) + '\n' )

	f.close()


#--------------------------------------------------------------------------------------------------------------------------------------
def parallel_calculate_normalisedEvents(kmer, fast5Files, poreModelFile, progress):
#	For hairpin training data - the large number of kmers makes it best to use parallel processing.  Uses a 5mer model to
#	calculate the shift and scale pore-specific parameters for each individual read in a list of fast5 files
#	ARGUMENTS
#       ---------
#	- kmer: redundant 6mer to identify the reads we're normalising
#	  type: string
#	- fast5Files: list of fast5 files whose events should be normalised
#	  type: list of strings
#	- poreModelFile: path to a pore model file.  This should be a 5mer model in the ONT format
#	  type: string
#	- progress: shows the kmer we're on to give an idea of progress
#	  type: tuple of length two
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

	sys.stdout.write("\rNormalising for shift and scale... " + str(progress[0]) + " of " + str(progress[1]))
	sys.stdout.flush()
	
	return (kmer, allNormalisedReads)


#--------------------------------------------------------------------------------------------------------------------------------------
def import_HairpinTrainingData(reference, BAMrecords, poreModelFile, redundant_A_Loc, position, readsThreshold):
#	Used to import training data from a hairpin primer of the form 5'-...NNNBNNN....NNNANNN...-3'.
#	Creates a map from kmer (string) to a list of lists, where each list is comprised of events from a read
#	First reads a BAM file to see which reads (readIDs, sequences) aligned to the references based on barcoding.  Then finds the fast5 files
#	that they came from, normalises the events to a pore model, and returns the list of normalised events.
#	ARGUMENTS
#       ---------
#	- reference: reference string from Osiris import_reference
#	  type: string
#	- bamFile: a BAM file from the alignment
#	  type: string
#	- poreModelFile: ONT model file for 5mers that can be used to normalised for shift and scale
#	  type: string
#	- redundant_A_Loc: location of the redundant A that is the reverse complement of BrdU (starting from 0)
#	  type: int
#	- position: where the analogue is in the 7mer (1&2, 3&4, 5&6)
#	  type: string
#	- readsThreshold: disregard a NNNANNN 7mer that only has a number of high quality reads below this threshold
#	  type: int
#	OUTPUTS
#       -------
#	- kmer2normalisedReads: a dictionary that takes a kmer string as a key and outputs a list of lists, where each list gives the normalised events from an individual read
#	  type: dictionary

	#count and print number of entries in bam file to stdout
	numRecords = len(BAMrecords)
	print str(numRecords) + ' records in BAM file.'

	if position == '12':
		startLB = redundant_A_Loc-5
		startUB = redundant_A_Loc-1
		endLB = redundant_A_Loc+6
		endUB = redundant_A_Loc+10
		Aloc = 1
	elif position == '34':
		startLB = redundant_A_Loc-7
		startUB = redundant_A_Loc-3
		endLB = redundant_A_Loc+4
		endUB = redundant_A_Loc+8
		Aloc = 3
	elif position == '56':
		startLB = redundant_A_Loc-9
		startUB = redundant_A_Loc-5
		endLB = redundant_A_Loc+2
		endUB = redundant_A_Loc+6
		Aloc = 5

	#build up the map that takes each indiviudal 7mer to a list of fast5 files that produced the reads
	kmer2Files = {}
	for record in BAMrecords:

		sequence = record.query_sequence
		readID = record.query_name

		#grab the part of the sequence that's flanked by start and end.  there may be more than one candidate.
		candidates = []
		start = reference[startLB:startUB] #four bases on the 5' end of the NNNANNN domain
		end = reference[endLB:endUB] #four bases on the 3' end of the NNNANNN domain
		start_indices = [s.start() for s in re.finditer('(?=' + start + ')', sequence)] #find all (possibly overlapping) indices of start using regular expressions
		end_indices = [s.start() for s in re.finditer('(?=' + end + ')', sequence)] #same for end
		for si in start_indices:
			si = si + len(start)
			for ei in end_indices:
				if ei > si:
					candidate = sequence[si:ei] #grab the subsequence between the start and end index
					if len(candidate) == 7 and candidate[Aloc] == 'A': #consider it a candidate if it's a 7mer and has an A in the middle
						candidates.append(candidate)

		#only add the read to the map if we're sure that we've found exactly one correct redundant 7mer, and its reverse complement is in the sequence
		if len(candidates) == 1:
			idx_brdu = sequence.find(reverseComplement(candidates[0]))
			idx_a = sequence.find(candidates[0])
			if idx_brdu != -1 and idx_brdu < idx_a:
				if candidates[0] in kmer2Files:
					kmer2Files[candidates[0]] += [readID]
				else:
					kmer2Files[candidates[0]] = [readID]

	#if a kmer has a number of associated reads that is below the minimum number of reads we need to train on, remove that kmer from the dictionary
	filteredKmer2Files = {}
	for key in kmer2Files:
		if len(kmer2Files[key]) >= readsThreshold:
			filteredKmer2Files[key] = kmer2Files[key]
	del kmer2Files

	#do the parallel processing, where a kmer (and associated reads) is given to each core.  use the maximum number of cores available
	normalisedReadsTuples = Parallel(n_jobs = multiprocessing.cpu_count())(delayed(parallel_calculate_normalisedEvents)(key, filteredKmer2Files[key], poreModelFile, (i,len(filteredKmer2Files))) for i, key in enumerate(filteredKmer2Files))

	#reshape the list of tuples from parallel processing into a dictionary
	kmer2normalisedReads = {}
	for entry in normalisedReadsTuples:
		kmer2normalisedReads[entry[0]] = entry[1]

	return kmer2normalisedReads


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
	os.system('bwa mem -t '+str(threads)+' -k4 -W4 -r2 -A1 -B1 -L3 -E0 '+pathToReference+' reads.fasta | samtools view -Sb - | samtools sort - alignments.sorted') 
	os.system('samtools index alignments.sorted.bam')

	sam_file = pysam.Samfile('alignments.sorted.bam')
	reference = import_reference(pathToReference)
	BAMrecords = []

	if position == '12':
		analogueLoc = reference.find('NTNNNNN')
		adenineLoc = reference.find('NNNNNAN')
		redundantALoc = adenineLoc + 6
	elif position == '34':
		analogueLoc = reference.find('NNNTNNN')
		adenineLoc = reference.find('NNNANNN')
		redundantALoc = adenineLoc + 4
	elif position == '56':
		analogueLoc = reference.find('NNNNNTN')
		adenineLoc = reference.find('NANNNNN')
		redundantALoc = adenineLoc + 2
	else:
		print 'Exiting with error.  Invalid argument passed to -p or --position.'
		splashHelp()

	for record in sam_file:

		if (record.aend is None) or (record.query_alignment_length is None) or record.is_reverse:
			continue

		#if the alignment indicates that this read doesn't have the critical regions, then ignore it and don't use it for training		
		if (record.reference_start > analogueLoc - 5) or (record.reference_end < adenineLoc + 10):
			continue
			
		BAMrecords.append(record)

	return (BAMrecords, redundantALoc)


#MAIN--------------------------------------------------------------------------------------------------------------------------------------
#parse arguments
args = sys.argv
a = parseArguments(args)

#do alignment QC
(BAMrecords, redundant_A_Loc) = alignAndSort(a.data, a.reference, a.position, a.threads)

#normalise the training data according to the ONT 5mer model
trainingData = import_HairpinTrainingData(import_reference(a.reference), BAMrecords, a.fiveMerModel, redundant_A_Loc, a.position, 40)

#write normalised reads to file in the .foh format
export_trainingDataToFoh( trainingData, a.outFoh )

