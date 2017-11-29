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
  -d,--data                 path to reads .fasta file or .bam file,
  -o,--output               path to the output detection .fdh file."""

	print s
	exit(0)


#--------------------------------------------------------------------------------------------------------------------------------------
def parseArguments(args):

	a = arguments()

	for i, argument in enumerate(args):
			
		if argument == '-d' or argument == '--data':
			a.reads = str(args[i+1])

		elif argument == '-o' or argument == '--output':
			a.outFdh = str(args[i+1])

		elif argument == '-h' or argument == '--help':
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
	raw_array = ( raw_array + offset ) * (rng/digitisation)

	f_hdf5.close()

	return raw_array


#--------------------------------------------------------------------------------------------------------------------------------------
def import_inVivoData(readsFile, outFile):

	f_out = open(outFile,'w')

	extension = readsFile.split('.')[-1:][0]

	#data is in fasta format
	if extension in ['fasta', 'fa']:
		f = open(readsFile,'r')
		g = f.readlines()
		f.close()

		for i, line in enumerate(g):

			#print progress
			displayProgress(i, len(g))

			#if you find a fasta header, it's either the first read or you've got a complete read with all bases
			if line[0] == '>':

				#if this isn't the first entry, write the filename, basecalls, and raw to the fdh file
				if i != 0:
					
					raw_array = readFast5Raw( filename )

					#write to fdh
					f_out.write('>' + filename + '\n' + bases + '\n')
					f_out.write( ' '.join(map(str,raw_array.tolist())) + '\n' )
				


				#regardless, set the new filename and empty the base string
				filename = line[1:].rstrip()
				bases = ''

			#otherwise, it's a basecall line (of which there may be several if we're sticking to 60char width)
			else:
				bases += line.rstrip()

		displayProgress(len(g), len(g))

	#data is in bam format from demultiplexing
	elif extension == 'bam':
		f = pysam.AlignmentFile(readsFile,'r')
		numOfRecords = f.count()
		counter = 0

		f = pysam.AlignmentFile(readsFile,'r')
		for record in f:
			filename = record.query_name
			raw_array = readFast5Raw( filename )
			bases = record.query_sequence

			#write to fdh
			f_out.write('>' + filename + '\n' + bases + '\n')
			f_out.write( ' '.join(map(str,raw_array.tolist())) + '\n' )
			
			displayProgress(counter, numOfRecords)
			counter += 1

		displayProgress(numOfRecords, numOfRecords)
		f.close()

	#bad or missing extension	
	else:
		print "Exiting with error: invalid extension in input file ", readsFile
		splashHelp()


	f_out.close()


#MAIN--------------------------------------------------------------------------------------------------------------------------------------
args = sys.argv
a = parseArguments(args)

import_inVivoData(a.reads, a.outFdh)

