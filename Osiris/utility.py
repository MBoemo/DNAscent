#!/usr/bin/env python


class BaseAnalogue():
#	class to house properties and data associated with a base analogue
#	PROPERTIES
#       ---------
#	- emissions: 
#	  type: dictionary (in the case of import_HairpinTrainingData) or list (in the case of import_FixedPosTrainingData)
#	- concentration: 
#	  type: positive float or int

	def __init__(self):
		self.emissions
		self.concentration = 0

	def set_emissions(self, trainingData):
		if trainingData == list or trainingData == dict:
			self.emissions = trainingData
		else:
			exit('Exiting: Training data must be a dictionary or a list.')

	def set_concentration(self, conc):
		if conc >= 0:
			self.concentration = conc
		else:
			exit('Exiting: Concentration must be a float greater than or equal to zero.')


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
