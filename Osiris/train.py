#!/usr/bin/env python

#----------------------------------------------------------
# Copyright 2017 University of Oxford
# Written by Michael A. Boemo (michael.boemo@path.ox.ac.uk)
#----------------------------------------------------------

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import sys

from data_IO import import_poreModel
from utility import reverseComplement, dynamicTimeWarping, warpPath, subsequenceDynamicTimeWarping
from build_model import build_TrainingHMM, build_softHMM
import numpy as np
from joblib import Parallel, delayed #for parallel processing
import multiprocessing #for parallel processing


def generateSignal(reference, poreModel):
#	generates a signal using the pore model.  this "ideal" signal is used in subsequence dynamic time warping
#	ARGUMENTS
#       ---------
#	- reference: output from import_reference from data_IO.py. 
#	  type: string
#	- poreModel: the imported ONT pore model (from import_poreModel)
#	  type: dictionary
#	OUTPUTS
#       -------
#	- signal: ideal signal produced by the reference
#	  type: list

	signal = []
	for i in range(0, len(reference)-5):
		signal.append(poreModel[reference[i:i+6]][0])
	return signal


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


def trainForFixedAnalogue(trainingData, reference, analoguePositions, poreModelFilename, threads):
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
#	- number of threads on which to run training
#	  type: int
#	OUTPUTS
#       -------
#	- analogueEmissions: keyed by a 6mer and returns a list of [trained_mean, trained_std]
#	  type: dictionary

	analogueEmissions = {}

	poreModel = import_poreModel(poreModelFilename)

	#build a training HMM based on the reference
	hmm = build_TrainingHMM(reference, poreModel)

	#train the HMM (Baum-Welch iterations) to the specified tolerance, using the specified number of threads	
	hmm.fit(trainingData, stop_threshold=1, max_iterations=50,verbose=True, n_jobs=threads)

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

	return (hmm, analogueEmissions)
			

def parallelSoftClip(key, Reads, reference, poreModel, progress, total):
#	this is a helper function for soft clipping in trainForContextAnalogue (below)
#	it does the soft clipping with dynamic time warping
#	ARGUMENTS
#       ---------
#	- key: BrdU 7mer from trainingData
#	  type: dictionary
#	- Reads: training reads for this 7mer
#	  type: list of lists (a keyed entry in trainingData)
#	- reference: reference sequence to build the HMM
#	  type: string
#	- poreModel: the imported ONT pore model (from import_poreModel)
#	  type: dictionary
#	- progress: 7mer that we're on (used to print progress)
#	  type: int
#	- total: number of 7mers in trainingData
#	  type: int
#	OUTPUTS
#       -------
#	- (key, truncatedReads): tuple of key and the corresponding training reads, which are now truncated by subsequence dynamic time warping
#	  type: tuple of (string, list of lists)

	#print progress
	sys.stdout.write("\rSoft clipping training data... " + str(progress) + " of " + str(total) )
	sys.stdout.flush()

	refLocal = reference
	revComp = reverseComplement(key)

	#replace the NNNANNN and NNNBNNN sequences in the reference with the 7mer for these reads
	refLocal = refLocal.replace('NNNANNN',key)
	refLocal = refLocal.replace('NNNTNNN',revComp)				

	#use the reference and pore model to generate an ideal signal that our region of interest should generate
	signalSnippet = generateSignal(refLocal, poreModel)

	#use subsequence dynamic time warping to truncate the reads so that they list events that are only in the region of interest (where NNNTNNN domain is)
	truncatedReads = []
	for i, read in enumerate(Reads):
		subseq = subsequenceDynamicTimeWarping( signalSnippet, read )
		if subseq is not None:
			truncatedReads.append(read[ max( [subseq[0], 0] ) : min( [subseq[1], len(read) - 1] ) ])

	return (key, truncatedReads)


def parallelTrainForSC(key, Reads, reference, poreModel, progress, total):
#	this is a helper function for soft clipping in trainForContextAnalogue (below)
#	ARGUMENTS
#       ---------
#	- key: BrdU 7mer from trainingData
#	  type: dictionary
#	- reads: training reads for this 7mer
#	  type: list of lists (a keyed entry in trainingData)
#	- reference: reference sequence to build the HMM
#	  type: string
#	- poreModel: the imported ONT pore model (from import_poreModel)
#	  type: dictionary
#	- progress: 7mer that we're on (used to print progress)
#	  type: int
#	- total: number of 7mers in trainingData
#	  type: int
#	OUTPUTS
#       -------
#	- analogueEmissions: keyed by a 6mer of the form NNBNNN or NNNBNN, and returns a list of [trained mean, trained standard deviation, number of training reads]
#	  type: dictionary

	#print progress
	sys.stdout.write("\rTraining for 6mer parameters... " + str(progress) + " of " + str(total) + ": ")
	sys.stdout.flush()
		
	refLocal = reference
	revComp = reverseComplement(key)

	#replace the NNNANNN and NNNBNNN sequences in the reference with the 7mer for these reads
	analogueLoc = reference.find('NNNTNNN')
	refLocal = refLocal.replace('NNNANNN',key)
	refLocal = refLocal.replace('NNNTNNN',revComp)

	#build a training HMM based on the reference
	hmm = build_softHMM(refLocal,poreModel)

	#train the HMM (Baum-Welch iterations) to the specified tolerance, using the specified number of threads	
	hmm.fit(Reads,edge_inertia=0,stop_threshold=1,n_jobs=1,verbose=False)

	analogueEmissions = {}

	#for each state
	for state in hmm.states:

		#ignore the states that have protected names (and do special things).  Note the addition of the garbage collection states, which we need for soft clipping
		if state.name not in ['None-end', 'None-start', 'startGarbageCollection', 'endGarbageCollection']: #these are the pomegranate protected names of start and end states
			
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
					
				analogueEmissions[kmer] = [ state.distribution.parameters[0], state.distribution.parameters[1], len(Reads) ] 

	return analogueEmissions


def trainForContextAnalogue(trainingData, reference, poreModelFilename, threads, softClip = False):
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
#	- threads: number of threads on which to run training
#	  type: int
#	- softClip: whether or not to use soft clipping to make the training run faster.  Soft clipping will only use a subsection of the read, around the ROI
#	  type: bool
#	OUTPUTS
#       -------
#	- analogueEmissions: keyed by a 6mer of the form NNBNNN or NNNBNN, and returns a list of [trained_mean, trained_std]
#	  type: dictionary

	
	#find analogue location in the reference
	analogueLoc = reference.find('NNNTNNN')
	poreModel = import_poreModel(poreModelFilename)

	#if we've specified to use soft clipping (finding the region of interest in the read with dynamic time warping)
	if softClip:

		#redefine reference as something shorter, and redefine analogueLoc as the location of the analogue in the truncated reference
		reference = reference[ max(analogueLoc-20, 0) : min(analogueLoc+26, len(reference)-1) ] 
		analogueLoc = reference.find('NNNTNNN')

		#do soft clipping (with dynamic time warping) in parallel
		results = Parallel(n_jobs=threads)(delayed(parallelSoftClip)(key, trainingData[key], reference, poreModel, progress, len(trainingData.keys())) for progress, key in enumerate(trainingData))

		#reshape results into dictionary, whereby training reads are keyed by 7mer (same as trainingData dictionary)
		truncatedTrainingData = {}
		for pair in results:
			truncatedTrainingData[pair[0]] = pair[1]

		#because we've truncated the read, we can get away with doing each 7mer in parallel
		#using the truncated reference and the truncated training data, do the training in parallel (one 7mer for each core)
		results = Parallel(n_jobs=threads)(delayed(parallelTrainForSC)(key, truncatedTrainingData[key], reference, poreModel, progress, len(truncatedTrainingData.keys())) for progress,key in enumerate(truncatedTrainingData))

		#reshape results into dictionary, whereby a BrdU-containing 6mer keys a vector of [mean, standard deviation, number of reads trained on]
		analogueEmissions = {}
		for emissDic in results:
			for key in emissDic:
				if key in analogueEmissions:
					if analogueEmissions[key][1] < emissDic[key][1]:
						pass
					else:
						analogueEmissions[key] = emissDic[key]
				else:
					analogueEmissions[key] = emissDic[key]

		return analogueEmissions

	#if soft clipping is not specified, then use the whole training read
	else:

		#this is the dictionary that we'll fill up with trained analogue emissions
		analogueEmissions = {}

		for i, key in enumerate(trainingData):

			#print progress
			sys.stdout.write("\rTraining for 6mer parameters... " + str(i) + " of " + str(len(trainingData.keys())) + ": ")
			sys.stdout.flush()
		
			#local reference for this loop, where we can replace the NNNTNNN and NNNANNN domains
			refLocal = reference
			revComp = reverseComplement(key)

			#replace the NNNANNN and NNNBNNN sequences in the reference with the 7mer for these reads
			refLocal = refLocal.replace('NNNANNN',key)
			refLocal = refLocal.replace('NNNTNNN',revComp)

			#build a training HMM based on the reference
			hmm = build_TrainingHMM(refLocal,poreModel)

			#train the HMM (Baum-Welch iterations) using the specified number of threads.  Here, we DON'T train transitions	
			hmm.fit(trainingData[key],edge_inertia=0,stop_threshold=0.1,n_jobs=threads)

			#go through the states of the trained HMM
			for state in hmm.states:

				#make sure the state is not protected.  these are the pomegranate protected names of start and end states
				if state.name not in ['None-end', 'None-start']:
 
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
						#if the kmer is already in the dictionary, check if the new one has a lower standard deviation than the one that's there
						#if it does, it's probably the better choice so use the new one instead.  Otherwise, stick with the old one.
						if kmer in analogueEmissions:
							if state.distribution.parameters[1] < analogueEmissions[kmer][1]:
								analogueEmissions[kmer] = [ state.distribution.parameters[0], state.distribution.parameters[1], len(trainingData[key]) ]	
						else:
							analogueEmissions[kmer] = [ state.distribution.parameters[0], state.distribution.parameters[1], len(trainingData[key]) ] 

		return analogueEmissions
