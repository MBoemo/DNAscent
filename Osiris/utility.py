#!/usr/bin/env python

#----------------------------------------------------------
# Copyright 2017 University of Oxford
# Written by Michael A. Boemo (michael.boemo@path.ox.ac.uk)
#----------------------------------------------------------


class BaseAnalogue():
#	class to house properties and data associated with a base analogue
#	PROPERTIES
#       ---------
#	- emissions: takes a 6mer and outputs emissions, output of trainForFixedAnalogue, trainForContextAnalogue, or import_poreModel (if you're importing an analogue pore model)
#	  type: dictionary
#	- concentration: 
#	  type: positive float or int

	def __init__(self, emis, conc):
		self.emissions = emis
		self.concentration = conc


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
			exit('Exiting: Illegal character in sequence.  Legal characters are A, T, G, and C.')

	#take the reverse of the complement and return it
	return revComp[::-1]
