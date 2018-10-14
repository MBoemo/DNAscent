#----------------------------------------------------------
# Copyright 2017 University of Oxford
# Written by Michael A. Boemo (michael.boemo@path.ox.ac.uk)
# This software is licensed under GPL-2.0.  You should have
# received a copy of the license with this software.  If
# not, please Email the author.
#----------------------------------------------------------


import sys
import matplotlib
matplotlib.use('Agg')
from matplotlib import pyplot as plt
import numpy as np
import matplotlib.mlab as mlab
from scipy import stats
import math

#outlier detection
from sklearn.cluster import DBSCAN
from sklearn import mixture

#--------------------------------------------------------------------------------------------------------------------------------------
class arguments:
	pass


#--------------------------------------------------------------------------------------------------------------------------------------
def splashHelp():
	s = """trainFromEventalign.py: Trains a BrdU mixture model from Nanopolish eventalign output.
To run trainFromEventalign.py, do:
  python trainFromEventalign.py [arguments]
Example:
  python trainFromEventalign.py -m /path/to/thymidine_model -e /path/to/nanopolish_eventalignment -o output_prefix
Required arguments are:
  -m,--ont_model            path to thymidine-only ONT model,
  -e,--alignment            path to Nanopolish eventalign output file,
  -o,--output               output prefix,
  -n,--maxReads             maximum number of reads to import from eventalign."""

	print s
	exit(0)


#--------------------------------------------------------------------------------------------------------------------------------------
def parseArguments(args):

	a = arguments()
	a.clipToMax = False
	a.maxReads = 1

	for i, argument in enumerate(args):
			
		if argument == '-m' or argument == '--ont_model':
			a.ont_model = str(args[i+1])

		elif argument == '-e' or argument == '--alignment':
			a.eventalign = str(args[i+1])

		elif argument == '-o' or argument == '--output':
			a.outFile = str(args[i+1])

		elif argument == '-n' or argument == '--maxReads':
			a.maxReads = int(args[i+1])
			a.clipToMax = True 

		elif argument == '-h' or argument == '--help':
			splashHelp()

	#check that required arguments are met
	if not hasattr( a, 'ont_model') or not hasattr( a, 'eventalign') or not hasattr( a, 'outFile') or not hasattr( a, 'maxReads'):
		splashHelp() 

	return a


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
def KLdivergence( mu1, sigma1, mu2, sigma2 ):

	return math.log(sigma2 / sigma1) + (sigma1**2 + (mu1 - mu2)**2)/(2.0*sigma2**2) - 0.5


#--------------------------------------------------------------------------------------------------------------------------------------
#parse arguments
args = sys.argv
a = parseArguments(args)

#eventalign output
sixmer2eventsBrdU = {}
f = open(a.eventalign,'r')
currentRead = ''
readCounter = 0
for line in f:

	splitLine = line.rstrip().split('\t')

	#ignore the header line
	if splitLine[0] == 'contig':
		continue

	readIndex = splitLine[3]
	if readIndex != currentRead:
		currentRead = readIndex
		readCounter += 1
		displayProgress(readCounter, a.maxReads)

	if a.clipToMax and readCounter == a.maxReads:
		break

	eventTime = float(splitLine[8])
	if eventTime < 0.002:
		continue

	sixmer = splitLine[9]
	if sixmer not in sixmer2eventsBrdU:
		sixmer2eventsBrdU[sixmer] = []
	elif len(sixmer2eventsBrdU[sixmer]) < 10000:
	#else:
		sixmer2eventsBrdU[sixmer].append( float(splitLine[6]) )
f.close()

#ONT pore model
model = {}
f = open(a.ont_model,'r')
for line in f:

	if line[0] == '#' or line[0] == 'k':
		continue
	
	splitLine = line.rstrip().split('\t')
	model[splitLine[0]] = [	float(splitLine[1]), float(splitLine[2]) ]
f.close()

f_out = open(a.outFile + '.brdu.model','w')
for i, key in enumerate(sixmer2eventsBrdU):

	if key == 'NNNNNN':
		continue
	
	if len( sixmer2eventsBrdU[key] ) > 200:
		x = np.linspace( np.mean(sixmer2eventsBrdU[key])-15, np.mean(sixmer2eventsBrdU[key])+15, 1000 )

		#noise reduction
		ar = np.array(sixmer2eventsBrdU[key])
		db = DBSCAN( min_samples= (0.025*len( sixmer2eventsBrdU[key] )) ).fit(ar.reshape(-1,1))
		outliers_filtered = []
		for j, label in enumerate(db.labels_):
			if label == -1:
				continue
			else:
				outliers_filtered.append(sixmer2eventsBrdU[key][j])

		#if we have enough events to train on after outlier detection
		if len(outliers_filtered) > 50:

			gmm = mixture.GMM(n_components=2, covariance_type='full')
			ar_filtered = np.array(outliers_filtered)
			out = gmm.fit(ar_filtered.reshape(-1,1))

			if abs(model[key][0] - out.means_[0]) < abs(model[key][0] - out.means_[1]):
				f_out.write( key + '\t' + str(model[key][0]) + '\t' + str(model[key][1]) + '\t' + str(out.weights_[1]) + '\t' + str(out.means_[1][0]) + '\t' + str(math.sqrt(out.covars_[1])) + '\t' + str(out.weights_[0]) + '\t' + str(out.means_[0][0]) + '\t' + str(math.sqrt(out.covars_[0])) + '\t' + str(KLdivergence( out.means_[1][0], math.sqrt(out.covars_[1]), model[key][0], model[key][1] ) ) + '\n')
			else:
				f_out.write( key + '\t' + str(model[key][0]) + '\t' + str(model[key][1]) + '\t' + str(out.weights_[0]) + '\t' + str(out.means_[0][0]) + '\t' + str(math.sqrt(out.covars_[0])) + '\t' + str(out.weights_[1]) + '\t' + str(out.means_[1][0]) + '\t' + str(math.sqrt(out.covars_[1])) + '\t' + str(KLdivergence( out.means_[0][0], math.sqrt(out.covars_[0]), model[key][0], model[key][1] ) ) + '\n')

f_out.close()


