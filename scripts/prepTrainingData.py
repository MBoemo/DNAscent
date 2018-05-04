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
def parallel_helper(sequence, bounds_r, bounds_q, path):

	#get the events from the fast5 file corresponding to this read, normalised to pA
	events = normaliseRead(path)

	return sequence + '\n' + bounds_r + '\n' + bounds_q + '\n' + events + '\n'


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

#first line in the output file should be the reference
f_out.write(reference + '\n')

numOfRecords = f_in.count()
f_out.write(str(numOfRecords) + '\n')

buffer_records = []
bufferMax = a.threads

#go through the BAM - when the buffer is full, get bounds and normalised events in parallel
f_in = pysam.Samfile(a.data,'r')
for counter, record in enumerate(f_in):

	#print progress to stdout
	displayProgress( counter, numOfRecords );
	
	#use the index to get the fast5 file path from the readID
	readRawPath =index[record.query_name]# record.query_name
	
	buffer_records.append( (record.query_sequence, str(record.reference_start) + ' ' + str(record.reference_end), str(record.query_alignment_start) + ' ' + str(record.query_alignment_end), readRawPath) )

	if len(buffer_records) > a.threads:

		out = Parallel(n_jobs=a.threads)(delayed(parallel_helper)(s,b_r,b_q,p) for s,b_r,b_q,p in buffer_records)

		for s in out:

			f_out.write(s)

		del buffer_records[:]
		gc.collect()

f_in.close()
f_out.close()
