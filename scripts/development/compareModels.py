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
import matplotlib
matplotlib.use('Agg')
from matplotlib import pyplot as plt
import numpy as np


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
			kmer2MeanStd[ splitLine[0] ] = [ float(splitLine[1]), float(splitLine[2]) ]
	g = None

	return kmer2MeanStd


#MAIN--------------------------------------------------------------------------------------------------------------------------------------
model1 = import_poreModel( sys.argv[1] )
model2 = import_poreModel( sys.argv[2] )

commonKeys = [key for key in model1 if key in model2]

diffs = []
for key in commonKeys:
	#if model1[key][1] < 3.5 and model2[key][1] < 3.5:
	diffs.append( model1[key][0] - model2[key][0] )

plt.hist(diffs,len(diffs)/10)
plt.xlabel('Difference Between Model Means (pA)')
plt.ylabel('Count')
plt.title('Comparison Between Models, N=' + str(len(diffs)))
plt.savefig('modelDifferencePlot.pdf')
print 'Mean: ',np.mean(diffs)
print 'Stdv: ',np.std(diffs)
