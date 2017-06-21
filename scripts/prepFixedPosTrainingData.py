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
  python prepFixedPosTrainingData.py -m /path/to/template_median68pA.5mer.model -d /path/to/alignment.bam -o output.foh -t 20
Required arguments are:
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
			
		elif argument == '-d' or argument == '--data':
			a.bamfile = str(args[i+1])

		elif argument == '-m' or argument == '--5mer-model':
			a.fiveMerModel = str(args[i+1])

		elif argument == '-o' or argument == '--output':
			a.outFoh = str(args[i+1])

		elif argument == '-h' or argument == '--help':
			splashHelp()

	return a


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

	for read in kmer2normalisedReads:
			
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


#MAIN--------------------------------------------------------------------------------------------------------------------------------------
#parse arguments
args = sys.argv
a = parseArguments(args)

#normalise the training data according to the ONT 5mer model
trainingData = import_FixedPosTrainingData(a.bamfile, a.fiveMerModel)

#write normalised reads to file in the .foh format
export_trainingDataToFoh( trainingData, a.outFoh )

