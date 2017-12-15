#----------------------------------------------------------
# Copyright 2017 University of Oxford
# Written by Michael A. Boemo (michael.boemo@path.ox.ac.uk)
# This software is licensed under GPL-3.0.  You should have
# received a copy of the license with this software.  If
# not, please Email the author.
#----------------------------------------------------------

#development script
#takes two Osiris models and plots a histogram of the difference between the mean for each 6mer
#how to run:
#	compareAtoB.py A_reference.fasta A_BAM.bam B_reference.fasta B_BAM.bam

import sys
import matplotlib
matplotlib.use('Agg')
from matplotlib import pyplot as plt
import numpy as np
import pysam


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
def import_AData(reference, BAMrecords, ROI, readsThreshold):
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

	analogueIndeces = range(ROI - 2, ROI + 9)

	#build up the map that takes each indiviudal 7mer to a list of fast5 files that produced the reads
	kmer2filename = {}

	#progress stuff
	numOfRecords = f.count()
	counter = 0
	fraction = 0.0
	failedROIQC = 0
	failedConcatenation = 0

	print "Binning BAM records..."

	f = pysam.Samfile(BAMrecords,'r')
	for record in f:

		#print progress
		counter += 1
		displayProgress( counter, numOfRecords)

		sequence = record.query_sequence

		if len(sequence) < int(1.1*len(reference)): #if this isn't a concatenated read...

			readID = record.query_name #read filename

			pairs = record.get_aligned_pairs(True,True) #tuples for each mapped position (pos-on-read,pos-on-ref,base-on-ref)

			analogueDomain = ['-']*len(analogueIndeces) 
			analogPosOnRead = [0]*len(analogueIndeces)

			#fill out analogueDomain and posOnRead using get_aligned_pairs() from pysam
			for p in pairs:
				if p[1] in analogueIndeces:
					analogueDomain[p[1] - (ROI - 2) ] = sequence[p[0]]
					analogPosOnRead[p[1] - (ROI - 2) ] = p[0]

			#make sure there aren't any indels
			indelFree = True
			for i,v in enumerate(analogPosOnRead[:-1]):
				if v != analogPosOnRead[i+1] - 1:
					indelFree = False

			brD = "".join(analogueDomain[2:9]).upper()
			LHS = "".join(analogueDomain[0:2]).upper()
			RHS = "".join(analogueDomain[9:]).upper()

			#if we've identified the analogue domain completely, keep this read to train on
			if ('-' not in analogueDomain) and (brD[3] == 'A') and (indelFree) and (LHS == reference[ROI-2:ROI]) and (RHS == reference[ROI+7:ROI+9]): 

				#append to dictionary
				if reverseComplement(brD) in kmer2filename:
					kmer2filename[reverseComplement(brD)] += [ ( readID, str(analogPosOnRead[0]) + ' ' + str(analogPosOnRead[-1]) ) ]
				else:
					kmer2filename[reverseComplement(brD)] = [ ( readID, str(analogPosOnRead[0]) + ' ' + str(analogPosOnRead[-1]) ) ]

			else:
				failedROIQC += 1
		else:
			failedConcatenation += 1

	displayProgress( numOfRecords, numOfRecords)
	f.close()
	return kmer2filename


#--------------------------------------------------------------------------------------------------------------------------------------
def import_BData(reference, BAMrecords, ROI, readsThreshold):
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

	analogueIndeces = range(ROI - 2, ROI + 9)

	#build up the map that takes each indiviudal 7mer to a list of fast5 files that produced the reads
	kmer2filename = {}

	#progress stuff
	numOfRecords = f.count()
	counter = 0
	fraction = 0.0
	failedROIQC = 0
	failedConcatenation = 0

	print "Binning BAM records..."

	f = pysam.Samfile(BAMrecords,'r')
	for record in f:

		#print progress
		counter += 1
		displayProgress( counter, numOfRecords)

		sequence = record.query_sequence

		if len(sequence) < int(1.1*len(reference)): #if this isn't a concatenated read...

			readID = record.query_name #read filename

			pairs = record.get_aligned_pairs(True,True) #tuples for each mapped position (pos-on-read,pos-on-ref,base-on-ref)

			analogueDomain = ['-']*len(analogueIndeces) 
			analogPosOnRead = [0]*len(analogueIndeces)

			#fill out analogueDomain and posOnRead using get_aligned_pairs() from pysam
			for p in pairs:
				if p[1] in analogueIndeces:
					analogueDomain[p[1] - (ROI - 2) ] = sequence[p[0]]
					analogPosOnRead[p[1] - (ROI - 2) ] = p[0]

			#make sure there aren't any indels
			indelFree = True
			for i,v in enumerate(analogPosOnRead[:-1]):
				if v != analogPosOnRead[i+1] - 1:
					indelFree = False

			brD = "".join(analogueDomain[2:9]).upper()
			LHS = "".join(analogueDomain[0:2]).upper()
			RHS = "".join(analogueDomain[9:]).upper()

			#if we've identified the analogue domain completely, keep this read to train on
			if ('-' not in analogueDomain) and (brD[3] == 'T') and (indelFree) and (LHS == reference[ROI-2:ROI]) and (RHS == reference[ROI+7:ROI+9]): 

				#append to dictionary
				if brD in kmer2filename:
					kmer2filename[brD] += [ ( readID, str(analogPosOnRead[0]) + ' ' + str(analogPosOnRead[-1]) ) ]
				else:
					kmer2filename[brD] = [ ( readID, str(analogPosOnRead[0]) + ' ' + str(analogPosOnRead[-1]) ) ]

			else:
				failedROIQC += 1
		else:
			failedConcatenation += 1

	displayProgress( numOfRecords, numOfRecords)
	f.close()
	return kmer2filename


#MAIN--------------------------------------------------------------------------------------------------------------------------------------
Areference = import_reference( sys.argv[1] )
Breference = import_reference( sys.argv[3] )

Amodel = import_AData( Areference, sys.argv[2], Areference.find('NNNANNN'), 40)
Bmodel = import_BData( Breference, sys.argv[4], Breference.find('NNNTNNN'), 40)

commonKeys = [key for key in Amodel if key in Bmodel]

x = []
y = []
for key in commonKeys:
	x.append( len(Amodel[key]) )
	y.append( len(Bmodel[key]) )

plt.scatter(x,y,alpha=0.3,edgecolors='none')
plt.ylim(ymin=0)
plt.xlim(xmin=0)
plt.xlabel('Number of training reads for Adenine')
plt.ylabel('Number of training reads for BrdU')
plt.axis('equal')
plt.title('Comparison Between BrdU and RC, N=' + str(len(x)))
plt.savefig('compareAtoB.pdf')
