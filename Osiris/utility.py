#!/usr/bin/env python

#----------------------------------------------------------
# Copyright 2017 University of Oxford
# Written by Michael A. Boemo (michael.boemo@path.ox.ac.uk)
#----------------------------------------------------------

import warnings
import math
import numpy as np


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


def hellingerDistance(mu1, std1, mu2, std2):
#	returns the Hellinger distance between normal distributions N(mu1,std1) and N(mu2,std2)

	h = 1 - math.sqrt( ( 2*std1*std2 )/( std1**2+std2**2 ) )*math.exp( ( -0.25*( mu1 - mu2 )**2 )/( std1**2 + std2**2 ) )
	h = math.sqrt( h )

	return h


def dynamicTimeWarping(x, y):
#	dyanmic time warping implementation for softclipping
#	ARGUMENTS
#       ---------
#	- x: vector of floats, which should be the generated desired events
#	  type: list of floats
#	- y: vector of floats, which should be all of the unfiltered events from the read
#	  type: list of floats
#	OUTPUTS
#       -------
#	- dtw: the dynamic time warping array for the two signals
#	  type: numpy array of floats
	
	#initialise matrix
	dtw = np.zeros( [len(x)+1, len(y)+1] )
	dtw[:,0] = np.inf
	dtw[0,:] = np.inf
	dtw[(0,0)] = 0
	
	#fill matrix at O(len(x)*len(y)) complexity
	for x_i, x_char in enumerate(x):
		x_i += 1
		for y_i, y_char in enumerate(y):
			y_i += 1
			dtw[ (x_i, y_i) ] = abs( x_char - y_char ) + min( [ dtw[ (x_i-1, y_i) ], dtw[ (x_i, y_i-1) ], dtw[ (x_i-1, y_i-1) ] ] )

	return dtw


def warpPath(dtw):
#	dyanmic time warping path for softclipping
#	ARGUMENTS
#       ---------
#	- x: "small" vector of floats, which should be the generated desired events
#	  type: list of floats
#	- y: "big" vector of floats, which should be all of the unfiltered events from the read
#	  type: list of floats
#	OUTPUTS
#       -------
#	- path: the best fit dynamic time warping path
#	  type: list of pairs of ints

	rows = dtw.shape[0] - 1
	cols = dtw.shape[1] - 1

	path = []

	i = rows
	j = cols

	while not ( i == 0 and j == 0):

		path.append( [i, j] )
		if i == 0:
			j -= 1
		elif j == 0:
			i -= 1
		else:
			m = np.argmin( [ dtw[(i-1,j-1)], dtw[(i-1,j)], dtw[ (i,j-1)] ] )
			if m == 0:
				j -= 1
				i -= 1
			elif m == 1:
				i -= 1
			elif m == 2:
				j -= 1
			else:
				exit(1)		

	return path


def subsequenceDynamicTimeWarping(x, y):
#	dyanmic time warping implementation for softclipping, specifically designed for subsequences
#	ARGUMENTS
#       ---------
#	- x: "small" vector of floats, which should be the generated desired events
#	  type: list of floats
#	- y: "big" vector of floats, which should be all of the unfiltered events from the read
#	  type: list of floats
#	OUTPUTS
#       -------
#	- (start,end): start index and end index of the portion of the long signal that best fits the short signal
#	  type: tuple of ints
	
	#initialise matrix, which is a little more complicated for this subsequence DTW variant
	dtw = np.zeros( [len(x), len(y)] )

	runningSum = 0
	for i in range(len(x)):
		runningSum += abs( x[i] - y[0] )
		dtw[i,0] = runningSum
	
	for j in range(len(y)):
		dtw[0,j] = abs( x[0] - y[j] )


	#fill matrix at O(len(x)*len(y)) complexity
	for x_i, x_char in enumerate(x[1:],start=1):
		for y_i, y_char in enumerate(y[1:],start=1):
			dtw[ (x_i, y_i) ] = abs( x_char - y_char ) + min( [ dtw[ (x_i-1, y_i) ], dtw[ (x_i, y_i-1) ], dtw[ (x_i-1, y_i-1) ] ] )

	#calculate path, starting from bStar
	bStar = np.argmin( dtw[ len(x) - 1, : ] )
	path = []
	i = len(x) - 1
	j = bStar

	#iterate back until one of i or j is 0
	while not ( i == 0 or j == 0 ):

		path.append( [i, j] )

		m = np.argmin( [ dtw[(i-1,j-1)], dtw[(i-1,j)], dtw[ (i,j-1)] ] )
		if m == 0:
			j -= 1
			i -= 1
		elif m == 1:
			i -= 1
		elif m == 2:
			j -= 1
		else:
			exit(1)		

	#take the start and end positions in the long signal and return them as a tuple
	start = path[-1][1]
	end = path[0][1]

	return (start,end)



