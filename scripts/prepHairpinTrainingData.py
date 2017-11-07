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
  -p,--position             position of analogue in training data (valid arguments are 1and2, 3and4, or 5and6),
  -n,--minimum              minimum threshold of reads that we're allowed to train on,
  -d,--data                 path to BAM file,
  -o,--output               path to the output training .foh or detection .fdh file."""

	print s
	exit(0)


#--------------------------------------------------------------------------------------------------------------------------------------
def parseArguments(args):

	a = arguments()

	for i, argument in enumerate(args):
			
		if argument == '-r' or argument == '--reference':
			a.reference = str(args[i+1])

		elif argument == '-d' or argument == '--data':
			a.data = str(args[i+1])

		elif argument == '-p' or argument == '--position':
			a.position = str(args[i+1])

		elif argument == '-o' or argument == '--output':
			a.outFoh = str(args[i+1])

		elif argument == '-n' or argument == '--minimum':
			a.minReads = int(args[i+1])

		elif argument == '-h' or argument == '--help':
			splashHelp()

	#check that required arguments are met
	if not hasattr( a, 'position') or not hasattr( a, 'data') or not hasattr( a, 'outFoh') or not hasattr( a, 'reference') or not hasattr( a, 'minReads'):
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
			
			rawStr = map( str, read[1] )
			f.write( read[0] + '\n' )
			f.write( ' '.join(rawStr) + '\n' )

	f.close()


#--------------------------------------------------------------------------------------------------------------------------------------
def import_HairpinTrainingData(reference, BAMrecords, ROI, readsThreshold):
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

	analogueIndeces = range(ROI[0] - 6, ROI[0] + 13)

	#build up the map that takes each indiviudal 7mer to a list of fast5 files that produced the reads
	kmer2raw = {}
	kmer2filename = {}

	#progress stuff
	numOfRecords = f.count()
	counter = 0
	fraction = 0.0

	f = pysam.Samfile(BAMrecords,'r')
	for record in f:

		#print progress
		counter += 1
		if round(float(counter)/float(numOfRecords),1) > fraction:
			sys.stdout.write('\rBinning BAM records: [' + '#'*int(fraction*10) + ' '*int(10-fraction*10) + '] ' + str(int(fraction*100)) + '%')
			sys.stdout.flush()
			fraction += 0.1

		sequence = record.query_sequence

		if len(sequence) < int(1.1*len(reference)): #if this isn't a concatenated read...

			readID = record.query_name #read filename

			pairs = record.get_aligned_pairs(True,True) #tuples for each mapped position (pos-on-read,pos-on-ref,base-on-ref)

			analogueDomain = ['-']*len(analogueIndeces) 
			analogPosOnRead = [0]*len(analogueIndeces)

			#fill out analogueDomain and posOnRead using get_aligned_pairs() from pysam
			for p in pairs:
				if p[1] in analogueIndeces:
					analogueDomain[p[1] - (ROI[0] - 6) ] = sequence[p[0]]
					analogPosOnRead[p[1] - (ROI[0] - 6) ] = p[0]

			#make sure there aren't any indels
			indelFree = True
			for i,v in enumerate(analogPosOnRead[:-1]):
				if v != analogPosOnRead[i+1] - 1:
					indelFree = False

			brD = "".join(analogueDomain[6:13]).upper()
			LHS = "".join(analogueDomain[0:6]).upper()
			RHS = "".join(analogueDomain[13:]).upper()

			#if we've identified the analogue domain completely, keep this read to train on
			if ('-' not in analogueDomain) and (analogueDomain[9] == 'T') and (indelFree) and (LHS == reference[ROI[0]-6:ROI[0]]) and (RHS == reference[ROI[0]+7:ROI[0]+13]): 

				#append to dictionary
				if brD in kmer2filename:
					kmer2filename[brD] += [readID] 
				else:
					kmer2filename[brD] = [readID]


	#only generate training data for kmers that have enough reads
	numOfRecords = len(kmer2filename.keys())
	fraction = 0.0
	for i, key in enumerate(kmer2filename):

		#print progress
		if round(float(i)/float(numOfRecords),1) > fraction:
			sys.stdout.write('\rConverting to pA: [' + '#'*int(fraction*10) + ' '*int(10-fraction*10) + '] ' + str(int(fraction*100)) + '%')
			sys.stdout.flush()
			fraction += 0.1


		if len(kmer2filename[key]) >= readsThreshold:

			kmer2raw[key] = []

			for filename in kmer2filename[key]:

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

				#append to dictionary
				kmer2raw[key] += [(sequence, raw_array.tolist())]

				f_hdf5.close()
	return kmer2raw


#MAIN--------------------------------------------------------------------------------------------------------------------------------------
#parse arguments
args = sys.argv
a = parseArguments(args)

reference = import_reference(a.reference)

if a.position == '1and2':
	analogueLoc = reference.find('NTNNNNN')
	adenineLoc = reference.find('NNNNNAN')

elif a.position == '3and4':
	analogueLoc = reference.find('NNNTNNN')
	adenineLoc = reference.find('NNNANNN')

elif a.position == '5and6':
	analogueLoc = reference.find('NNNNNTN')
	adenineLoc = reference.find('NANNNNN')

else:
	print 'Exiting with error.  Invalid argument passed to -p or --position.'
	splashHelp()

#normalise the training data according to the ONT 5mer model
trainingData = import_HairpinTrainingData(import_reference(a.reference), a.data, [analogueLoc, adenineLoc], a.minReads)

#write normalised reads to file in the .foh format
export_trainingDataToFoh( trainingData, a.outFoh )

print '\n'

