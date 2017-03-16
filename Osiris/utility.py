#!/usr/bin/env python

#----------------------------------------------------------
# Copyright 2017 University of Oxford
# Written by Michael A. Boemo (michael.boemo@path.ox.ac.uk)
#----------------------------------------------------------

import warnings


class BaseAnalogue():
#	class to house properties and data associated with a base analogue
#	PROPERTIES
#       ---------
#	- emissions: takes a 6mer and outputs emissions, output of trainForFixedAnalogue, trainForContextAnalogue, or import_poreModel (if you're importing an analogue pore model)
#	  type: dictionary
#	- concentration: 
#	  type: positive float or int

	def __init__(self, emis, conc):
		if type(emis) == dict:
			self.emissions = emis
		else:
			warnings.warn('Warning: First argument to BaseAnalogue constructor must be a dictionary.', Warning)
		if type(conc) == float and conc >= 0.0 and conc <= 1.0:
			self.concentration = conc
		else:
			warnings.warn('Warning: Second argument to BaseAnalogue constructor must be a positive float between 0.0 and 1.0.', Warning)


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
