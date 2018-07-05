#----------------------------------------------------------
# Copyright 2017 University of Oxford
# Written by Michael A. Boemo (michael.boemo@path.ox.ac.uk)
# This software is licensed under GPL-3.0.  You should have
# received a copy of the license with this software.  If
# not, please Email the author.
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
from joblib import Parallel, delayed
import multiprocessing

#--------------------------------------------------------------------------------------------------------------------------------------
class arguments:
	pass


#--------------------------------------------------------------------------------------------------------------------------------------
def splashHelp():
	s = """prepTrainingData.py: Osiris preprocessing script that will format ONT reads in the Osiris training/detection format.
To run prepTrainingData.py, do:
  python prepTrainingData.py [arguments]
Example:
  python prepTrainingData.py -r /path/to/reference.fasta -d /path/to/aligned.bam -o output.foh -t 20
Required arguments are:
  -r,--reference            path to the reference file,
  -i,--index                path to the index file made by index.py,
  -d,--data                 path to BAM file,
  -o,--output               path to the output training .foh or detection .fdh file.
Optional arguments are:
  -t,--threads              number of threads (default is 1 thread)."""

	print s
	exit(0)


#--------------------------------------------------------------------------------------------------------------------------------------
def parseArguments(args):

	a = arguments()
	a.threads = 1

	for i, argument in enumerate(args):
			
		if argument == '-r' or argument == '--reference':
			a.reference = str(args[i+1])

		elif argument == '-d' or argument == '--data':
			a.data = str(args[i+1])

		elif argument == '-i' or argument == '--index':
			a.indexfile = str(args[i+1])

		elif argument == '-o' or argument == '--output':
			a.outFoh = str(args[i+1])

		elif argument == '-t' or argument == '--threads':
			a.threads = int(args[i+1])

		elif argument == '-h' or argument == '--help':
			splashHelp()

	#check that required arguments are met
	if not hasattr( a, 'data') or not hasattr( a, 'outFoh') or not hasattr( a, 'reference') or not hasattr( a, 'indexfile'):
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

	f = open(filename,'r')
	g = f.readlines()
	f.close()

	reference = {}
	for line in g:

		if line[0] == '>':
			currentChromosome = line[1:].rstrip()
			reference[ currentChromosome ] = ''

		if line[0] != '>':
			reference[ currentChromosome ] += line.rstrip().upper()

	return reference


#--------------------------------------------------------------------------------------------------------------------------------------
def normaliseRead(filename):

	f_hdf5 = h5py.File(filename,'r')
	
	#get raw
	path = '/Raw/Reads'
	read_number = f_hdf5[path].keys()[0]
	raw_data = f_hdf5[path + '/' + read_number + '/Signal']
	raw_array = np.zeros(raw_data.len(),dtype='int16')
	raw_data.read_direct(raw_array)

	#get read-specific parameters to normalise from raw to pA
	scaling = f_hdf5['/UniqueGlobalKey/channel_id']
	digitisation = scaling.attrs.get('digitisation')
	offset = scaling.attrs.get('offset')
	rng = scaling.attrs.get('range')
	sample_rate = scaling.attrs.get('sampling_rate')

	#normalise to pA
	raw_array = ( raw_array + offset ) * (rng/digitisation)

	f_hdf5.close()

	return ' '.join(map(str,raw_array.tolist()))


#--------------------------------------------------------------------------------------------------------------------------------------
def parallel_helper(refSequence, basecall, path, pairs):

	#get the events from the fast5 file corresponding to this read, normalised to pA
	events = normaliseRead(path)

	return refSequence + '\n' + basecall + '\n' + events + '\n' + pairs + '\n'

#--------------------------------------------------------------------------------------------------------------------------------------
def fix_pairs( pairs ):

	fixed = []
	currentR = pairs[0][1]
	currentQ = pairs[0][0]
	fixed.append( (currentQ, currentR) )

	for p in pairs[1:]:

		if p[1] == currentR + 1:
			currentR = p[1]
			currentQ = p[0]
			fixed.append( (currentQ, currentR) )
		else:
			rollingR = p[1]
			currentQ = p[0]
			while rollingR != currentR:
				fixed.append( (currentQ, rollingR) )
				rollingR -= 1
			currentR = p[1]
			
	return fixed
				
#MAIN--------------------------------------------------------------------------------------------------------------------------------------

#parse arguments
args = sys.argv
a = parseArguments(args)

#load the index
f_index = open(a.indexfile,'r')
index = {}
for line in f_index:
	splitLine = line.rstrip().split('\t')
	index[splitLine[0]] = splitLine[1]
f_index.close()

#open input/output files
f_in = pysam.Samfile(a.data,'r')
f_out = open(a.outFoh,'w')

#load the reference
reference = import_reference(a.reference)

numOfRecords = 1000#f_in.count()
f_out.write(str(numOfRecords) + '\n')

buffer_records = []
bufferMax = a.threads

#go through the BAM - when the buffer is full, get bounds and normalised events in parallel
f_in = pysam.Samfile(a.data,'r')
for counter, record in enumerate(f_in):

	if record.is_reverse == False:
		counter += 1

		#print progress to stdout
		displayProgress( counter, numOfRecords );
	
		#use the index to get the fast5 file path from the readID
		readRawPath =index[record.query_name]

		#get the subsequence of the reference that this read mapped to
		refAlignedSequence = reference[ record.reference_name ][ record.reference_start: record.reference_end ]

		#get the Albacore basecall
		basecall = record.query_sequence

		#aligned pairs
		pairs = fix_pairs(record.get_aligned_pairs(True,False))
		pairs_str = ''
		correctiveFactor = pairs[0][1]
		if correctiveFactor != record.reference_start:
			print "WARNING"


		for p in pairs:
			pairs_str += str(p[0]) + ' ' + str(p[1] - correctiveFactor) + ' '
	
		buffer_records.append( (refAlignedSequence, basecall, readRawPath, pairs_str) )

		if len(buffer_records) > a.threads:

			out = Parallel(n_jobs=a.threads)(delayed(parallel_helper)(s,b,p,pairs) for s,b,p,pairs in buffer_records)

			for s in out:

				f_out.write(s)

			del buffer_records[:]
			gc.collect()

		if counter == 1000:
			break

f_in.close()
f_out.close()
