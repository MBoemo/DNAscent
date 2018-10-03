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


#------------------------------------------------------------------------------------------------------------------------------------------
def divergence(mu1,sig1,mu2,sig2):

	return math.log(sig2/sig1) + ( math.pow(sig1,2) + math.pow(mu1-mu2,2))/(2*math.pow(sig2,2)) - 0.5


#------------------------------------------------------------------------------------------------------------------------------------------
def import_poreModel(filename,c1,c2):
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
			kmer2MeanStd[ splitLine[0] ] = ( float(splitLine[c1]), float(splitLine[c2]), float(splitLine[1]), float(splitLine[2]) )
	g = None

	return kmer2MeanStd


#MAIN--------------------------------------------------------------------------------------------------------------------------------------
plt.figure()

model = import_poreModel( sys.argv[1], 4, 5 )
diffs = []
for key in model:
	mu1, sig1, mu2, sig2 = model[key]
	diffs.append( divergence(mu1,sig1,mu2,sig2) )

plt.subplot(2,1,1)
plt.hist(diffs,50,linewidth=0)
plt.tick_params(axis='both', which='major', labelsize=8)
plt.tick_params(axis='both', which='minor', labelsize=8)
plt.ylabel('Count')
plt.title('Expected Log Likelihood of BrdU to ONT Model')
axes = plt.gca()
axes.set_xlim([0,3.5])

model = import_poreModel( sys.argv[1], 7, 8 )
diffs = []
for key in model:
	mu1, sig1, mu2, sig2 = model[key]
	diffs.append( divergence(mu1,sig1,mu2,sig2) )

plt.subplot(2,1,2)
plt.hist(diffs,50,linewidth=0)
plt.tick_params(axis='both', which='major', labelsize=8)
plt.tick_params(axis='both', which='minor', labelsize=8)
plt.xlabel('Expected Log Likelihood')
plt.ylabel('Count')
plt.title('Expected Log Likelihood of Trained Thymidine to ONT Model')
axes = plt.gca()
axes.set_xlim([0,3.5])


plt.tight_layout()
plt.savefig('modelDifferencePlot.pdf')
