#----------------------------------------------------------
# Copyright 2017 University of Oxford
# Written by Michael A. Boemo (michael.boemo@path.ox.ac.uk)
# This software is licensed under GPL-3.0.  You should have
# received a copy of the license with this software.  If
# not, please Email the author.
#----------------------------------------------------------

#for hairpin training data
#preps data for C++ Osiris from reads run on a 1D flow cell, basecalled using 1D settings


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
	s = """prepHairpinTrainingData.py: Osiris preprocessing script that will format ONT reads in the Osiris training/detection format.
To run prepHairpinTrainingData.py, do:
  python prepHairpinTrainingData.py [arguments]
Example:
  python prepHairpinTrainingData.py -r /path/to/reference.fasta -p 34 -m /path/to/template_median68pA.5mer.model -d /path/to/reads -o output.foh -t 20
Required arguments are:
  -r,--reference            path to reference file in fasta format,
  -p,--position             position of analogue in training data (valid arguments are 23, 34, or 45),
  -d,--data                 path to BAM file,
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

	#check that required arguments are met
	if not hasattr( a, 'fiveMerModel') or not hasattr( a, 'position') or not hasattr( a, 'data') or not hasattr( a, 'outFoh') or not hasattr( a, 'reference'):
		splashHelp() 

	return a


#--------------------------------------------------------------------------------------------------------------------------------------
def reverseComplement(sequence):
#	takes a DNA sequence and returns its reverse complement
#	ARGUMENTS
#       ---------
#	- sequence: DNA sequence
#	  type: string
#	OUTPUTS
#       -------
#	- revComp: reverse complement DNA sequence
#	  type: string

	#seq to upper case
	sequence = sequence.upper()

	#build the complement
	revComp = ''
	for char in sequence:
		if char == 'A':
			revComp += 'T'
		elif char == 'T':
			revComp += 'A'
		elif char == 'C':
			revComp += 'G'
		elif char == 'G':
			revComp += 'C'
		#guard against illegal characters
		else:
			warnings.warn('Warning: Illegal character in sequence.  Legal characters are A, T, G, and C.', Warning)

	#take the reverse of the complement and return it
	return revComp[::-1]


#--------------------------------------------------------------------------------------------------------------------------------------
def import_poreModel(filename):
#	takes the filename of an ONT pore model file and returns a map from kmer (string) to [mean,std] (list of floats)
#	ARGUMENTS
#       ---------
#	- filename: path to an ONT model file
#	  type: string
#	OUTPUTS
#       -------
#	- kmer2MeanStd: a map, keyed by a kmer, that returns the model mean and standard deviation signal for that kmer
#	  type: dictionary

	f = open(filename,'r')
	g = f.readlines()
	f.close()

	kmer2MeanStd = {}
	for line in g:
		if line[0] != '#' and line[0:4] != 'kmer': #ignore the header
			splitLine = line.split('\t')
			kmer2MeanStd[ splitLine[0] ] = [ float(splitLine[1]), float(splitLine[2]) ]
	g = None

	return kmer2MeanStd


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

		A = np.zeros((3,3))
		b = np.zeros((3,1))
		for event in Events:
			 if float(event[7]) > 0.3: #if there's a high probability (>30%) that the 5mer model called by Metrichor was the correct one
				model_5mer = event[4]
				event_mean = float(event[0])
				event_time = float(event[2])
				model_mean = kmer2MeanStd[model_5mer][0]
				model_std = kmer2MeanStd[model_5mer][1]
				
				#update matrix A
				A[0,0] += 1/(model_std**2)
				A[0,1] += 1/(model_std**2)*model_mean
				A[0,2] += 1/(model_std**2)*event_time
				A[1,1] += 1/(model_std**2)*model_mean**2
				A[1,2] += 1/(model_std**2)*model_mean*event_time
				A[2,2] += 1/(model_std**2)*event_time**2

				#update vector b
				b[0] += 1/(model_std**2)*event_mean
				b[1] += 1/(model_std**2)*event_mean*model_mean
				b[2] += 1/(model_std**2)*event_mean*event_time

		#use symmetry of A
		A[1,0] = A[0,1]
		A[2,0] = A[0,2]
		A[2,1] = A[1,2]

		#solve Ax = b to find shift and scale
		x = np.linalg.solve(A,b)
		shift = x[0][0]
		scale = x[1][0]
		drift = x[2][0]

		#go through the same events as before and normalise them to the pore model using scale and shift
		normalisedEvents = []
		for event in Events:
			if float(event[7]) > 0.30: #if there's a high probability (>30%) that the 5mer model called by Metrichor was the correct one
				event_mean = float(event[0])
				normalisedEvents.append( (event_mean - shift - drift*event_time)/scale )

		allNormalisedReads.append(normalisedEvents)

		f.close()

	sys.stdout.write("\rNormalising for shift and scale... " + str(progress[0]) + " of " + str(progress[1]))
	sys.stdout.flush()
	
	return (kmer, allNormalisedReads)


#--------------------------------------------------------------------------------------------------------------------------------------
def import_HairpinTrainingData(reference, BAMrecords, poreModelFile, ROI, readsThreshold):
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
	
	f = pysam.Samfile(BAMrecords,'r')

	analogueIndeces = range(ROI[0], ROI[0] + 7)
	adenineIndeces = range(ROI[1], ROI[1] + 7)

	#build up the map that takes each indiviudal 7mer to a list of fast5 files that produced the reads
	kmer2Files = {}
	for record in f:

		sequence = record.query_sequence

		if len(sequence) < int(1.1*len(reference)): #if this isn't a concatenated read...

			readID = record.query_name #read filename

			pairs = record.get_aligned_pairs(True,True) #tuples for each mapped position (pos-on-read,pos-on-ref,base-on-ref)

			adenineDomain = ['-']*7 #fill this up as we identify bases for the 7mer

			for p in pairs:
				if p[1] in adenineIndeces: #if the reference position is in the redundant adenine domain
					adenineDomain[p[1] - ROI[1]] = sequence[p[0]] #add the base from the read that mapped to it to the adenine 7mer


			adD = "".join(adenineDomain).upper() #make a string out of the adenineDomain list
			if ('-' not in adenineDomain) and (adD[3] == 'A'): #if we've identified the adenine domain completely

				if adD in kmer2Files:
					kmer2Files[adD] += [readID]
				else:
					kmer2Files[adD] = [readID]

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


#MAIN--------------------------------------------------------------------------------------------------------------------------------------
#parse arguments
args = sys.argv
a = parseArguments(args)

reference = import_reference(a.reference)

if a.position == '23':
	analogueLoc = reference.find('NNTNNNN')
	adenineLoc = reference.find('NNNNANN')

elif a.position == '34':
	analogueLoc = reference.find('NNNTNNN')
	adenineLoc = reference.find('NNNANNN')

elif a.position == '45':
	analogueLoc = reference.find('NNNNTNN')
	adenineLoc = reference.find('NNANNNN')

else:
	print 'Exiting with error.  Invalid argument passed to -p or --position.'
	splashHelp()

#normalise the training data according to the ONT 5mer model
trainingData = import_HairpinTrainingData(import_reference(a.reference), a.data, a.fiveMerModel, [analogueLoc, adenineLoc], 40)

#write normalised reads to file in the .foh format
export_trainingDataToFoh( trainingData, a.outFoh )

print '\n'

