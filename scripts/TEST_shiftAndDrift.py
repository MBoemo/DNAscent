import numpy as np
import sys
import warnings
import h5py
import pysam
import re
import os
import gc
import math


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
def calculate_normalisedEvents(fast5Files):


	#open the 5mer model and make a map from the 5mer (string) to [mean,std] (list)
	kmer2MeanStd = import_poreModel('/home/nieduszynski/michael/development/Osiris/pore_models/template_median68pA.5mer.model')

	#now iterate through all the relevant fast5 files so we only need to open the model file once
	allNormalisedReads = []
	for f5File in fast5Files:

		#use 5mer model to calculate shift, scale, drift, and var to normalise events for the pore
		f = h5py.File(f5File,'r')
		path = '/Analyses/Basecall_1D_000/BaseCalled_template/Events'
		Events = f[path]

		A = np.zeros((2,2))
		b = np.zeros((2,1))
		for event in Events:
			 if float(event[7]) > 0.3: #if there's a high probability (>30%) that the 5mer model called by Metrichor was the correct one
				model_5mer = event[4]
				event_mean = float(event[0])
				event_time = float(event[2])
				model_mean = kmer2MeanStd[model_5mer][0]
				model_std = kmer2MeanStd[model_5mer][1]

				#update matrix A
				A[0,0] += 1/(model_std**2)
				A[0,1] += 1/(model_std**2)*model_mean
				A[1,1] += 1/(model_std**2)*model_mean**2

				#update vector b
				b[0] += 1/(model_std**2)*event_mean
				b[1] += 1/(model_std**2)*event_mean*model_mean

		#use symmetry of A
		A[1,0] = A[0,1]


		#solve Ax = b to find shift and scale
		x = np.linalg.solve(A,b)
		shift = x[0][0]
		scale = x[1][0]

		print "shift: ",shift
		print "scale: ", scale


		#go through the same events as before and normalise them to the pore model using scale and shift
		normalisedEvents = []
		for event in Events:
			if float(event[7]) > 0.30: #if there's a high probability (>30%) that the 5mer model called by Metrichor was the correct one
				event_mean = float(event[0])
				#ne = (event_mean - shift - drift*event_time)/scale
				ne = (event_mean - shift )/scale
				if ne > 50 and ne < 130:
					normalisedEvents.append( ne )


		allNormalisedReads.append(normalisedEvents)

		f.close()
#--------------------------------------------------------------------------------------------------------------------------------------
f = pysam.Samfile( sys.argv[1], 'r')

filenames = []

for record in f:
	
	filenames.append( record.query_name )

calculate_normalisedEvents(filenames)
