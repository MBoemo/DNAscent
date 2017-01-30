#!/usr/bin/env python

#----------------------------------------------------------
# Copyright 2017 University of Oxford
# Written by Michael A. Boemo (michael.boemo@path.ox.ac.uk)
#----------------------------------------------------------

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

from data_IO import import_poreModel
from utility import reverseComplement
from build_model import build_TrainingHMM
import numpy as np


def callAnalogue(hmm, readEvents):
#	For plotting the difference between the trained mean/std and the pore model mean/std
#	ARGUMENTS
#       ---------
#	- hmm: Hidden Markov Model object, built by build_RandIncHMM from build_model.py
#	  type: pomegranate object
#	- readEvents: normalised events from an Oxford Nanopore reads
#	  type: list
#	OUTPUTS
#       -------
#	- calledAnaloguePositions: list of positions on the reference where the model has detected a base analogue
#	  type: list

	calledAnaloguePositions = []

	#run the viterbi algorithm on the readEvents.  the best path of states is given by vpath and quality is given by logp
	logp, vpath = hmm.viterbi(readEvents)

	for entry in vpath:
		if entry[1].name != 'None-end' and entry[1].name != 'None-start': #protected pomegranate names for start and end states
			splitName = entry[1].name.split('_') #state info: [branch (T or B), state type (I, D, M1, etc.), 'pos', position on reference]
			
			#if the state is a BrdU state and it was a match (either M1 or M2) then call a base analogue at that position
			if splitName[0] == 'B' and splitName[1] in ['M1','M2']:
				calledAnaloguePositions.append(int(splitName[3]))

	#in theory, there could be both M1 and M2 matches for a single HMM module, so remove duplicates
	calledAnaloguePositions = list(set(calledAnaloguePositions))

	return calledAnaloguePositions


def trainingImprovements(hmm, reference, modelFile, kmerName):
#	For plotting the difference between the trained mean/std and the pore model mean/std
#	ARGUMENTS
#       ---------
#	- hmm: the trained hidden markov model from trainForFixedAnalogue or trainForContextAnalogue
#	  type: dictionary
#	- reference: output from import_reference from data_IO.py. 
#	  type: string
#	- modelFile: an ONT pore model filename for a 6mer model
#	  type: string
#	- kmerName: mainly for trainForContextAnalogue, the name of the context kmer (just for plotting and filenames)
#	  type: string
#	OUTPUTS
#       -------
#	- meanImprovement, stdImprovement: lists of the difference between the trained mean or std and model mean or std at each position on the reference
#	  type: lists
	
	#make a map from kmer to model emissions
	modelEmissions = import_poreModel(modelFile)

	#make a map reference position to kmer
	posToKmer = {}
	for i, char in enumerate(reference[:-6]):
		posToKmer[i] = reference[i:i+6]

	meanImprovement = [0]*(i+1)
	stdImprovement = [0]*(i+1)
	for state in hmm.states:
		if state.name != 'None-end' and state.name != 'None-start':
			stateInfo = state.name.split('_') #state info: [branch (T or B), state type (I, D, M1, etc.), 'pos', position on reference]
			if stateInfo[1] in ['M1','M2']:
				kmer = posToKmer[int(stateInfo[3])]
				meanImprovement[int(stateInfo[3])] = abs(state.distribution.parameters[0] - modelEmissions[kmer][0])
				stdImprovement[int(stateInfo[3])] = abs(state.distribution.parameters[1] - modelEmissions[kmer][1])

	#plot mean improvement
	x = range(len(meanImprovement))
	plt.figure()
	plt.bar(x,meanImprovement,width=1.0)
	plt.title('Training Improvement for '+kmerName)
	plt.ylabel('|Model Mean - Trained Mean| (pA)')
	plt.xlabel('Position on Reference (bp)')
	plt.savefig('trainingImprovement_'+kmerName+'_mu.pdf')
	plt.close()

	#plot standard deviation improvement
	plt.figure()
	plt.bar(x,stdImprovement,width=1.0)
	plt.ylabel('|Model Std - Trained Std| (pA)')
	plt.xlabel('Position on Reference (bp)')
	plt.savefig('trainingImprovement_'+kmerName+'_std.pdf')
	plt.close()

	return meanImprovement, stdImprovement


def trainForFixedAnalogue(trainingData, reference, analoguePositions, poreModelFilename, tolerance, threads):
#	For training data with one base analogue in a fixed context - no N's
#	Creates a map from kmer (string) to a list that has the kmer's mean and standard deviation
#	ARGUMENTS
#       ---------
#	- trainingData: data, should be the output of import_FixedPosTrainingData from data_IO.py
#	  type: list
#	- reference: output from import_reference from data_IO.py
#	  type: string
#	- AnaloguePositions: position of the analogues in the reference, where indexing starts from 0 (first base is given index 0). If there is only one analogue, should be a list of one element, e.g. [55]
#	  type: list
#	- poreModelFilename: an ONT pore model filename for a 6mer model
#	  type: string
#	- toerance for convergence of Baum-Welch iterations
#	  type: float
#	- number of threads on which to run training
#	  type: int
#	OUTPUTS
#       -------
#	- analogueEmissions: keyed by a 6mer and returns a list of [trained_mean, trained_std]
#	  type: dictionary

	analogueEmissions = {}

	#build a training HMM based on the reference
	hmm = build_TrainingHMM(reference,poreModelFilename)

	#train the HMM (Baum-Welch iterations) to the specified tolerance, using the specified number of threads	
	hmm.fit(trainingData,stop_threshold=tolerance,n_jobs=threads)

	for state in hmm.states:
		if state.name != 'None-end' and state.name != 'None-start': #these are the pomegranate protected names of start and end states
			stateInfo = state.name.split('_') #state info: [branch (T or B), state type (I, D, M1, etc.), 'pos', position on reference]
			stateLoc = int(stateInfo[3])
			for analogueLoc in analoguePositions:
				if (stateInfo[1] in ['M1','M2']) and (stateLoc in range(analogueLoc-5, analogueLoc+1)):
					kmer = reference[stateLoc:stateLoc+6]				
				
					#replace the T in the reference with a B for the base analogue
					kmer = kmer[0:analogueLoc-stateLoc] + 'B' + kmer[analogueLoc-stateLoc+1:]
					
					#dictionary, keyed by the analogue 6mer, that returns the trained mean and trained standard deviation for the 6mer
					analogueEmissions[kmer] = [ state.distribution.parameters[0], state.distribution.parameters[1] ] 

	return analogueEmissions
			

def trainForContextAnalogue(trainingData, reference, poreModelFilename, tolerance, threads):
#	Creates a map from kmer (string) to list of the kmer's mean and standard deviation 
#	First reads a BAM file to see which reads (readIDs, sequences) aligned to the references based on barcoding.  Then finds the fast5 files
#	that they came from, normalises the events to a pore model, and returns the list of normalised events.
#	ARGUMENTS
#       ---------
#	- trainingData: data, should be the output of import_HairpinTrainingData from data_IO.py
#	  type: dictionary
#	- reference: output from import_reference from data_IO.py.  Should have NNNTNNN for BrdU 7mer and NNNANNN for redundant 7mer
#	  type: string
#	- poreModelFilename: an ONT pore model filename for a 6mer model
#	  type: string
#	- toerance for convergence of Baum-Welch iterations
#	  type: float
#	- number of threads on which to run training
#	  type: int
#	OUTPUTS
#       -------
#	- analogueEmissions: keyed by a 6mer of the form NNBNNN or NNNBNN, and returns a list of [trained_mean, trained_std]
#	  type: dictionary

	analogueEmissions = {}
	for key in trainingData:

		refLocal = reference
		revComp = reverseComplement(key)

		#calculate analogue location in the reference
		analogueLoc = refLocal.find('NNNTNNN')

		#replace the NNNANNN and NNNBNNN sequences in the reference with the 7mer for these reads
		refLocal = refLocal.replace('NNNANNN',key)
		refLocal = refLocal.replace('NNNTNNN',revComp)

		#build a training HMM based on the reference
		hmm = build_TrainingHMM(refLocal,poreModelFilename)

		#train the HMM (Baum-Welch iterations) to the specified tolerance, using the specified number of threads	
		hmm.fit(trainingData[key],stop_threshold=tolerance,n_jobs=threads)

		for state in hmm.states:
			if state.name != 'None-end' and state.name != 'None-start': #these are the pomegranate protected names of start and end states
				stateInfo = state.name.split('_') #state info: [branch (T or B), state type (I, D, M1, etc.), 'pos', position on reference]
				stateLoc = int(stateInfo[3])
				if (stateInfo[1] in ['M1','M2']) and (stateLoc in [analogueLoc, analogueLoc + 1]):
					#replace the T in the reference with a B for the base analogue
					kmer = refLocal[stateLoc:stateLoc+6]
					if stateLoc == analogueLoc:
						kmer = kmer[0:3] + 'B' + kmer[4:]
					elif stateLoc == analogueLoc + 1:
						kmer = kmer[0:2] + 'B' + kmer[3:]
					else:
						exit('Exiting: Problem in trainForAnalogue in identifying analogue location.')
					
					#dictionary, keyed by the analogue 6mer, that returns the trained mean and trained standard deviation for the 6mer
					analogueEmissions[kmer] = [ state.distribution.parameters[0], state.distribution.parameters[1] ] 
		
		#PLOTTING: this is to plot training improvements per kmer as a pdf , if you want to (uncomment if this is desired)
		#NOTE: this produces two plots per kmer, so depending on what you're doing, it might produce a lot of plots
		#row_mu, row_std = trainingImprovements(hmm, refLocal, poreModelFilename, key)
		
		#PLOTTING: this block of code creates numpy arrays that can be plotted if you want to visualise a heatmap of mean, std improvements across a number of kmers
		#if first:
		#	A = np.array(row_mu)
		#	B = np.array(row_std)
		#	first = False
		#else:
		#	A = np.vstack((A,row_mu))
		#	B = np.vstack((B,row_std))

	return analogueEmissions
