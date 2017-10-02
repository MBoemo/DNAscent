#----------------------------------------------------------
# Copyright 2017 University of Oxford
# Written by Michael A. Boemo (michael.boemo@path.ox.ac.uk)
#----------------------------------------------------------

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
	s = """prepDetectionData.py: Osiris preprocessing script that will format ONT reads in the Osiris detection format.
To run prepDetectionData.py, do:
  python prepDetectionData.py [arguments]
Example:
  python prepDetectionData.py -m /path/to/template_median68pA.5mer.model -d /path/to/reads.fasta -o output.fdh -t 20
Required arguments are:
  -d,--data                 path to reads .fasta file or .bam file,
  -m,--5mer-model           path to 5mer pore model file (provided by ONT) to normalise reads,
  -o,--output               path to the output detection .fdh file.
Optional arguments are:
  -t,--threads              number of threads (default is 1 thread)."""

	print s
	exit(0)


#--------------------------------------------------------------------------------------------------------------------------------------
def parseArguments(args):

	a = arguments()

	#defaults
	a.threads = 1

	for i, argument in enumerate(args):
		if argument == '-t' or argument == '--threads':
			a.threads = int(args[i+1])
			
		elif argument == '-d' or argument == '--data':
			a.reads = str(args[i+1])

		elif argument == '-m' or argument == '--5mer-model':
			a.fiveMerModel = str(args[i+1])

		elif argument == '-o' or argument == '--output':
			a.outFdh = str(args[i+1])

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
def normaliseSingleRead(fast5File, bases, kmer2MeanStd, progress, total):
#	For fixed position training data - small enough to be done in serial.  Uses a 5mer model to calculate the shift and scale
#	pore-specific parameters for each individual read in a list of fast5 files
#	ARGUMENTS
#       ---------
#	- fast5File: single fast5 file
#	  type: string
#	- poreModelFile: path to a pore model file.  This should be a 5mer model in the ONT format
#	  type: string
#	OUTPUTS
#       -------
#	- allNormalisedRead: normalised events for the read
#	  type: list

	sys.stdout.write("\rNormalising in vivo reads... " + str(progress) + " of " + str(total))
	sys.stdout.flush()

	#use 5mer model to calculate shift, scale, drift, and var to normalise events for the pore
	f = h5py.File(fast5File,'r')
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
			event_time = float(event[2])
			normalisedEvent = (event_mean - shift - drift*event_time)/scale
			if normalisedEvent > 0.0:
				normalisedEvents.append(normalisedEvent)

	f.close()
	
	return (fast5File, bases, normalisedEvents)


#--------------------------------------------------------------------------------------------------------------------------------------
def import_inVivoData(readsFile, fiveMerPoreModelFile, n_threads):


	extension = readsFile.split('.')[-1:][0]
	linesDic = {}

	if extension in ['fasta', 'fa']:
		f = open(readsFile,'r')
		g = f.readlines()
		f.close()

		for i, line in enumerate(g):
			if line[0] == '>':
				if i != 0:
					linesDic[filename] = bases				

				filename = line[1:].rstrip()
				bases = ''
			else:
				bases += line.rstrip()

	elif extension == 'bam':
		f = pysam.AlignmentFile(readsFile,'r')

		for record in f:
			linesDic[record.query_name] = record.query_sequence


	poreModel = import_poreModel(fiveMerPoreModelFile)

	results = Parallel(n_jobs = n_threads)(delayed(normaliseSingleRead)(filename, linesDic[filename], poreModel, i, len(linesDic.keys())) for i, filename in enumerate(linesDic))
	
	return results


#--------------------------------------------------------------------------------------------------------------------------------------
def export_toFdh( fileTuples, filename ):
#	takes reads from a run, aligns them to a reference, and separates the resulting bam file by each reference
#	ARGUMENTS
#       ---------
#	- fileTuples: result of import_inVivoData which is a list of tuples of the form (fast5 filename, basecalls, normalised events)
#	  type: list of tuples
#	- filename: output .fdh filename
#	  type: string
#	OUTPUTS
#       -------
#	- exports read observation data in the .fdh format, to be read by C++ Osiris
	
	f = open(filename,'w')

	for threeTuple in fileTuples:

		#if any of the fields are empty, skip this read
		if ('' in [threeTuple[0], threeTuple[1]]) or (threeTuple[2] == []):
			continue
		
		f.write('>' + threeTuple[0] + '\n' + threeTuple[1] + '\n')
			
		eventsStr = map(str, threeTuple[2])

		f.write( ' '.join(eventsStr) + '\n' )

	f.close()


#MAIN--------------------------------------------------------------------------------------------------------------------------------------
#parse arguments
args = sys.argv
a = parseArguments(args)

results = import_inVivoData(a.reads, a.fiveMerModel, a.threads)

export_toFdh( results, a.outFdh )

