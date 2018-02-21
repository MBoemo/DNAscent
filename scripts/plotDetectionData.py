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
import numpy as np
import matplotlib
matplotlib.use('Agg')
from matplotlib import pyplot as plt

#--------------------------------------------------------------------------------------------------------------------------------------
class arguments:
	pass


#--------------------------------------------------------------------------------------------------------------------------------------
def splashHelp():
	s = """plotDetectionData.py: Osiris postprocessing script that will plot analogue calls against reference position.
To run plotDetectionData.py, do:
  python plotDetectionData.py [arguments]
Example:
  python plotDetectionData.py -d /path/to/detectionFile.detect -b /path/to/alignment.bam
Required arguments are:
  -d,--detection            path to detection file produced by Osiris detect,
  -b,--bam                  path to bam file used to make the .fdh detection input file,
  -r,--reference            path to reference fasta file."""

	print s
	exit(0)


#--------------------------------------------------------------------------------------------------------------------------------------
def parseArguments(args):

	a = arguments()
	a.threads = 1

	for i, argument in enumerate(args):

		if argument == '-d' or argument == '--detection':
			a.detect = str(args[i+1])

		elif argument == '-b' or argument == '--bam':
			a.bam = str(args[i+1])

		elif argument == '-r' or argument == '--reference':
			a.reference = str(args[i+1])


	#check that required arguments are met
	if not hasattr( a, 'detect') or not hasattr( a, 'bam') or not hasattr( a, 'bam'):
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


#MAIN----------------------------------------------------------------------------------------------------------------------------------
args = sys.argv
a = parseArguments(args)

#import reference 
reference = import_reference(a.reference)

#make a map between filename and aligned pairs
f_bam = pysam.Samfile(a.bam,'r')
fname2align = {}
for record in f_bam:
	
	fname2align[record.query_name] = record.get_aligned_pairs()
f_bam.close()

#use the alignment to adjust the 
f_detect = open(a.detect,'r')
g = f_detect.readlines()
allCalls = []
callVector = [-1]*len(reference)
first = True
BrdUcalls = 0.0
ThymidineCalls = 0.0

for line in g:
	if line[0] == '>':
		pairs = fname2align[line[1:].rstrip()]

		if not first:

			allCalls.append(callVector)
			callVector = [-1]*len(reference)

		first = False
	else:
		#get the position on the reference
		splitLine = line.rstrip().split('\t')
		callFloat = float(splitLine[1])
		positionOnRead = int(splitLine[0])
		positionOnRef = None
		for p in pairs:
			if p[0] == positionOnRead:
				positionOnRef = p[1]
				break
		
		if positionOnRef is not None:
			if callFloat >= 2.5:
				callVector[positionOnRef] = 1
				BrdUcalls += 1.0
			else:
				callVector[positionOnRef] = 0
				ThymidineCalls += 1.0			

na = np.array(allCalls)

fig = plt.figure()
cax = plt.matshow(na,cmap='hot')
fig.colorbar(cax)

plt.savefig('out.pdf')

print "Percent BrdU calls: ", str(100.0*BrdUcalls/(BrdUcalls+ThymidineCalls))

