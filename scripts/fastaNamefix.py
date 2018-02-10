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
from joblib import Parallel, delayed
import multiprocessing

#--------------------------------------------------------------------------------------------------------------------------------------
class arguments:
	pass


#--------------------------------------------------------------------------------------------------------------------------------------
def splashHelp():
	s = """fastaNamefix.py: Osiris preprocessing script that create reads.fasta with read names as fast5 paths.
To run fastaNamefix.py, do:
  python fastaNamefix.py [arguments]
Example:
  python fastaNamefix.py -f /path/to/fastq -i path/to/index.osiris
Required arguments are:
  -f,--fastq                path to fastq from basecalling,
  -i,--index                path to Osiris index file."""

	print s
	exit(0)


#--------------------------------------------------------------------------------------------------------------------------------------
def parseArguments(args):

	a = arguments()
	a.threads = 1

	for i, argument in enumerate(args):
		if argument == '-f' or argument == '--fastq':
			a.fastq = str(args[i+1])

		elif argument == '-i' or argument == '--index':
			a.index = str(args[i+1])

		elif argument == '-h' or argument == '--help':
			splashHelp()
		elif argument[0] == '-':
			splashHelp()

	#check that required arguments are met
	if not hasattr( a, 'fastq') or not hasattr( a, 'index'):
		splashHelp() 

	return a


#MAIN--------------------------------------------------------------------------------------------------------------------------------------
args = sys.argv
a = parseArguments(args)

#load the index
f_index = open('index.osiris','r')
index = {}
for line in f_index:
	splitLine = line.rstrip().split('\t')
	index[splitLine[0]] = splitLine[1]
f_index.close()

#go through the fastq
f_fastq = open(a.fastq,'r')
f_out = open('reads.fasta','w')
grapLine = False
for line in f_fastq:
	if line[0] == '@':
		name = line.rstrip().split(' ')[0]
		name = name[1:]
		grabLine = True
	elif grabLine:
		seq = line.rstrip()
		grabLine = False
		f_out.write('>' + index[name] + '\n' + seq + '\n')
f_fastq.close()
f_out.close()
