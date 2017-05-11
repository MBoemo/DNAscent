#!/usr/bin/env python

#----------------------------------------------------------
# Copyright 2017 University of Oxford
# Written by Michael A. Boemo (michael.boemo@path.ox.ac.uk)
#----------------------------------------------------------

from build_model import build_RandIncHMM
from utility import reverseComplement
from data_IO import import_poreModel
from joblib import Parallel, delayed # for parallel processing
import multiprocessing
from train import generateSignal
from utility import hellingerDistance, subsequenceDynamicTimeWarping
import numpy as np
import math


def callHairpin(kmer2normalisedReads, reference, poreModelFilename, analogueObj):
#	Test function for calling analogues in hairpin primer data (data of the type that we used for training)
#	ARGUMENTS
#       ---------
#	- kmer2normalisedReads: output of import_HairpinTrainingData from data_IO.py
#	  type: dictionary
#	- reference: output of import_reference from data_IO.py.  Has NNNTNNN and NNNANNN domains.
#	  type: string
#	- poreModelFilename: filename for the ONT 6mer poremodel over {A,T,G,C}
#	  type: string
#	- analogueObj: Osiris analogue object from utility.py
#	  type: Osiris object
#	OUTPUTS
#       -------
#	- calledAnaloguePositions: a list of tuples, where the first entry is the NNNANNN kmer and the second entry 
#	  is a list of lists, where each entry is a list of called analogue positions in a single read
#	  type: list

	#filter the analogue emissions so that it only has the emissions that we can distinguish from thymidine
	print "Kmers in analogue: " + str(len(analogueObj.emissions.keys()))
	ontModel = import_poreModel(poreModelFilename)
	filteredEmissions = {}
	for key in analogueObj.emissions:
		if analogueObj.emissions[key][1] < 3 and hellingerDistance( analogueObj.emissions[key][0], analogueObj.emissions[key][1], ontModel[key.replace('B','T')][0], ontModel[key.replace('B','T')][1] ) > 0.2:

			filteredEmissions[key] = analogueObj.emissions[key]

	analogueObj.emissions = filteredEmissions
	print "Filtered to distinguishable Kmers: " + str(len(analogueObj.emissions.keys()))
	
	kmerAndCalls = []
	filtered = {}
	references = {}

	for key in kmer2normalisedReads:
		
		refLocal = reference
		revComp = reverseComplement(key)

		#calculate analogue location in the reference
		analogueLoc = refLocal.find('NNNTNNN')

		#replace the NNNANNN and NNNBNNN sequences in the reference with the 7mer for these reads
		refLocal = refLocal.replace('NNNANNN',key)
		refLocal = refLocal.replace('NNNTNNN',revComp)

		#if we have analogue emissions for the 7mer in this hairpin, try to call the analogue
		candidiate = revComp[0:3] + 'B' + revComp[4:]

		if (candidiate[0:6] in analogueObj.emissions) and (candidiate[1:] in analogueObj.emissions):

			references[key] = refLocal
			filtered[key] = kmer2normalisedReads[key]

			#uncomment for serial
			#calledAnaloguePositions = softCallAnalogue(kmer2normalisedReads[key], refLocal, poreModelFilename, analogueObj)

			#uncomment for serial
			#kmerAndCalls.append( (candidiate, calledAnaloguePositions) )

	#comment for serial
	results = Parallel( n_jobs = multiprocessing.cpu_count() )( delayed( softCallAnalogue )( filtered[ key ], references[ key ], poreModelFilename, analogueObj ) for key in filtered )

	return results


def hardCallAnalogue( normalisedReads, reference, poreModelFilename, analogueObj ):
#	Calls the positions of a base analogue in a single read
#	ARGUMENTS
#       ---------
#	- normalisedReads: list of normalised events for a single read
#	  type: list
#	- reference: output of import_reference from data_IO.py.
#	  type: string
#	- poreModelFilename: filename for the ONT 6mer poremodel over {A,T,G,C}
#	  type: string
#	- analogueObj: Osiris analogue object from utility.py
#	  type: Osiris object
#	OUTPUTS
#       -------
#	- calledAnaloguePositions: list of positions on the reference where the model has detected a base analogue
#	  type: list

	calledAnaloguePositions = []

	#build an HMM that can detect the random incorporation of a base analogue that replaces thymidine
	hmm = build_RandIncHMM(reference, poreModelFilename, analogueObj)

	for read in normalisedReads:

		#run the viterbi algorithm on the normalisedReads.  the best path of states is given by vpath and quality is given by logp
		logp, vpath = hmm.viterbi(read)

		analoguesInThisRead = []
		for entry in vpath:
			if entry[1].name != 'None-end' and entry[1].name != 'None-start': #protected pomegranate names for start and end states
				splitName = entry[1].name.split('_') #state info: [branch (T or B), state type (I, D, M1, etc.), 'pos', position on reference]
			
				#if the state is a BrdU state and it was a match (either M1 or M2) then call a base analogue at that position
				if splitName[0] == 'B' and splitName[1] in ['M1','M2']:
					analoguesInThisRead.append(int(splitName[3]))

		#in theory, there could be both M1 and M2 matches for a single HMM module, so remove duplicates
		analoguesInThisRead = list(set(analoguesInThisRead))

		calledAnaloguePositions.append( (logp, analoguesInThisRead ) )

	return calledAnaloguePositions

def logSum(x,y):
	if math.isnan(x) or math.isnan(y):
		if math.isnan(x) and math.isnan(y):
			return float('nan')
		elif math.isnan(x):
			return y
		else:
			return x
	else:
		if x > y:
			return x + np.log( 1 + np.exp(y - x))
		else:
			return y + np.log(1 + np.exp(x - y))
	

def logProd(x,y):
	if math.isnan(x) or math.isnan(y):
		return float('nan')
	else:
		return x + y	


def softCallAnalogue( normalisedReads, reference, poreModelFilename, analogueObj ):
#	Calls the positions of a base analogue in a single read
#	ARGUMENTS
#       ---------
#	- normalisedReads: list of normalised events for a single read
#	  type: list
#	- reference: output of import_reference from data_IO.py.
#	  type: string
#	- poreModelFilename: filename for the ONT 6mer poremodel over {A,T,G,C}
#	  type: string
#	- analogueObj: Osiris analogue object from utility.py
#	  type: Osiris object
#	OUTPUTS
#       -------
#	- calledAnaloguePositions: list of positions on the reference where the model has detected a base analogue
#	  type: list

	calledAnaloguePositions = []

	#build an HMM that can detect the random incorporation of a base analogue that replaces thymidine
	detectHMM = build_RandIncHMM(reference, poreModelFilename, analogueObj)
	analoguePositions = detectHMM[0]
	hmm = detectHMM[1]
	stateIndices = { state:i for i, state in enumerate( hmm.states ) }
	
	for read in normalisedReads:

		probabilityAtPos = {}
		Tprob = {}
		Bprob = {}

		#initialise the T and B probabilities at each position to 0
		for position in analoguePositions:
			Tprob[position] = float('nan')
			Bprob[position] = float('nan')

		#run the forward algorithm on this read
		#forwardLattice = hmm.forward(read)
		#backwardLattice = hmm.backward(read)
		forwardLattice = hmm.predict_log_proba(read)

		#for each state in the HMM...
		for state in hmm.states:
			
			if state.name != 'None-end' and state.name != 'None-start': #protected pomegranate names for start and end states
				stateSplit = state.name.split('_')
				posOnRef = int(stateSplit[3])
				stateName = stateSplit[1]
				TorB = stateSplit[0]

				#if this state is at a position on the reference where we can make a BrdU call
				if (posOnRef in analoguePositions) and (stateName in ['M1','M2']):
					rowInForward = stateIndices[state]
					
					if TorB == 'T':
						for i,j in enumerate(forwardLattice[:,rowInForward]):
							Tprob[posOnRef] = logSum(Tprob[posOnRef], j)#logProd(j,backwardLattice[i,rowInForward]))

					elif TorB == 'B':
						for i,j in enumerate(forwardLattice[:,rowInForward]):
							Bprob[posOnRef] = logSum(Bprob[posOnRef], j)#logProd(j,backwardLattice[i,rowInForward]))

		for key in Tprob:

			probabilityAtPos[key] = Bprob[key] - Tprob[key]

		calledAnaloguePositions.append(probabilityAtPos)

	return calledAnaloguePositions
