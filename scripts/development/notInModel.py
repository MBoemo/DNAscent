#----------------------------------------------------------
# Copyright 2017 University of Oxford
# Written by Michael A. Boemo (michael.boemo@path.ox.ac.uk)
# This software is licensed under GPL-3.0.  You should have
# received a copy of the license with this software.  If
# not, please Email the author.
#----------------------------------------------------------

#development script
#takes two Osiris models and plots a histogram of the difference between the mean for each 6mer

import sys
import math
import matplotlib
matplotlib.use('Agg')
from matplotlib import pyplot as plt
import numpy as np
from itertools import product

#------------------------------------------------------------------------------------------------------------------------------------------
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
			kmer2MeanStd[ splitLine[0] ] = ( float(splitLine[1]), float(splitLine[2]), float(splitLine[3]), float(splitLine[4]) )
	g = None

	return kmer2MeanStd


#MAIN--------------------------------------------------------------------------------------------------------------------------------------
model = import_poreModel( sys.argv[1] )
twoMers = [''.join(i) for i in product(['A','T','G','C'],repeat=2)]
threeMers = [''.join(i) for i in product(['A','T','G','C'],repeat=3)]

BrdUpos3_4 = []
for lhs in twoMers:
	for rhs in threeMers:
		BrdUpos3_4.append( lhs+'B'+rhs )

for lhs in threeMers:
	for rhs in twoMers:
		BrdUpos3_4.append( lhs+'B'+rhs )	

for item in BrdUpos3_4:
	if item not in model:
		print item
	


