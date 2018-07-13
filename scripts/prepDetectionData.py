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
  python prepDetectionData.py -d /path/to/reads.fasta -o output.fdh
Required arguments are:
  -r,--reference            path to the reference file,
  -i,--index                path to the index file made by index.py,
  -d,--data                 path to alignment bam file,
  -o,--output               path to the output detection .fdh file."""

	print s
	exit(0)


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
def parseArguments(args):

	a = arguments()

	for i, argument in enumerate(args):
			
		if argument == '-d' or argument == '--data':
			a.reads = str(args[i+1])

		elif argument == '-r' or argument == '--reference':
			a.reference = str(args[i+1])

		elif argument == '-i' or argument == '--index':
			a.index = str(args[i+1])

		elif argument == '-o' or argument == '--output':
			a.outFdh = str(args[i+1])

		elif argument == '-h' or argument == '--help':
			splashHelp()

	#check that required arguments are met
	if not hasattr( a, 'reads') or not hasattr( a, 'outFdh') or not hasattr( a, 'reference') or not hasattr( a, 'index'):
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
def readFast5Raw( filename ):

	#read fast5 file
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
	raw_array = (raw_array + offset)*(rng / digitisation)

	f_hdf5.close()

	return raw_array


#--------------------------------------------------------------------------------------------------------------------------------------
def import_inVivoData(readsFile, outFile):

	f_out = open(outFile,'w')

	#data is in bam format from demultiplexing
	f = pysam.AlignmentFile(readsFile,'r')
	numOfRecords = f.count()
	counter = 0

	#write the foh header
	f_out.write( str(numOfRecords) + '\n' )

	f = pysam.AlignmentFile(readsFile,'r')
	for record in f:
		filename = record.query_name
		raw_array = readFast5Raw( filename )
		bases = record.query_sequence

		#write to fdh
		f_out.write( '>' + filename + '\n' + bases + '\n' )
		f_out.write( ' '.join(map(str,raw_array.tolist())) + '\n<\n' )
			
		displayProgress(counter, numOfRecords)
		counter += 1

	displayProgress(numOfRecords, numOfRecords)
	f.close()
	f_out.close()


#MAIN--------------------------------------------------------------------------------------------------------------------------------------
args = sys.argv
a = parseArguments(args)

#load the index
f_index = open(a.indexfile,'r')
index = {}
for line in f_index:
	splitLine = line.rstrip().split('\t')
	index[splitLine[0]] = splitLine[1]
f_index.close()

#load the reference
reference = import_reference(a.reference)



import_inVivoData(a.reads, a.outFdh)

