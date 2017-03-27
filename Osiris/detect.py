#!/usr/bin/env python

#----------------------------------------------------------
# Copyright 2017 University of Oxford
# Written by Michael A. Boemo (michael.boemo@path.ox.ac.uk)
#----------------------------------------------------------

from build_model import build_RandIncHMM
from utility import reverseComplement


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

	kmerAndCalls = []

	for key in kmer2normalisedReads:
		
		refLocal = reference
		revComp = reverseComplement(key)

		#calculate analogue location in the reference
		analogueLoc = refLocal.find('NNNTNNN')

		#replace the NNNANNN and NNNBNNN sequences in the reference with the 7mer for these reads
		refLocal = refLocal.replace('NNNANNN',key)
		refLocal = refLocal.replace('NNNTNNN',revComp)

		calledAnaloguePositions = callAnalogue(kmer2normalisedReads[key], refLocal, poreModelFilename, analogueObj)

		kmerAndCalls.append( (key, calledAnaloguePositions) )

	return kmerAndCalls


def callAnalogue( normalisedReads, reference, poreModelFilename, analogueObj ):
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

		calledAnaloguePositions.append(analoguesInThisRead)

	return calledAnaloguePositions
