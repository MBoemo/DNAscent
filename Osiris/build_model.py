#!/usr/bin/env python

#----------------------------------------------------------
# Copyright 2017 University of Oxford
# Written by Michael A. Boemo (michael.boemo@path.ox.ac.uk)
#----------------------------------------------------------

from pomegranate import *
import numpy as np
import json
import time
import warnings
from data_IO import import_reference, import_poreModel


def build_RandIncHMM(refSequence, poreModel, analogue):
#	Builds a HMM from a reference sequence with a topology appropriate for detecting randomly incorporated base analogue.
#	ARGUMENTS
#       ---------
#	- refSequence: reference sequence from import_reference
#	  type: string
#	- poreModel: the imported ONT pore model (from import_poreModel)
#	  type: dictionary
#	- analogue: an Osiris object created by utility.py that contains the concentration and emission data for a base analogue
#	  type: Osiris object
#	OUTPUTS
#       -------
#	- hmm: a hidden Markov model with forks for base analogues
#	  type: pomegranate object
	
	#positions where we create an analogue module
	analoguePositions = []

	#set the hmm class
	hmm = HiddenMarkovModel() 

	#get length of reference
	refLength = len(refSequence)

	#add analogue emissions to the emissions for A-T-G-C 6mers from the ONT pore model file
	poreModel.update(analogue.emissions)
	analogueConc = analogue.concentration

	#two-dimensional array for the states, where the columns are the positions in the reference
	thymidineStates = [[0 for x in range(refLength)] for y in range(6)] 

	#Dictionary for analogue states.  Keyed by an integer position on the reference, where an analogue is about to come into the 4th position.
	#It returns a list lists containing pomegranate state objects.
	#Example: analogueStates[i] = [ [analogue states for position i+1], [analogue states for position i+2] ]
	analogueStates = {} 

	#Initial transitions within modules (internal transitions)
	internalSS2M1 = 0.97
	internalSS2M2 = 0.03
	internalD2I = 0.01
	internalI2I = 0.50
	internalI2SS = 0.49
	internalM12M1 = 0.51
	internalM12SE = 0.49
	internalM22M2 = 0.97
	internalM22SE = 0.03
	internalSE2I = 0.01

	#Initial transitions between modules (external transitions)
	externalD2D = 0.85
	externalD2SS = 0.14
	externalI2SS = 0.01
	externalSE2D = 0.12
	externalSE2SS = 0.87

	#this is really just a safety thing... get the last index in the iterator
	for last_i,x in enumerate(refSequence[:-5]): 
		pass
		
	###########################################################################################################################
	# Add States to Model
	
	#Create the HMM states.  Iterate through the reference sequence, and make the repeating HMM module for each position in the sequence.
	#Stop 7 characters before the end, because we'll need to reach forward by 7 to handle the branching.
	for i, char in enumerate(refSequence[:-5]): 

		###################################################
		# Handle Thymidine Branch States

		#grab the appropriate emission probabilities for the 6mer that we're on
		presentEmissions = poreModel[refSequence[i:i+6]]

		thymidineStates[0][i] = State( None, name='T_SS_pos_'+str(i) )                                                          #0 SS
		thymidineStates[1][i] = State( None, name='T_D_pos_'+str(i) )                                                           #1 D
		thymidineStates[2][i] = State( UniformDistribution(30, 130), name='T_I_pos_'+str(i) )                                   #2 I  - uniformly distributed on 30pA to 130pA
		thymidineStates[3][i] = State( None, name='T_M1_pos_'+str(i) )                                                          #3 M1 - we'll tie this to M2 in a minute
		thymidineStates[4][i] = State( NormalDistribution(presentEmissions[0], presentEmissions[1]), name='T_M2_pos_'+str(i) )  #4 M2 
		thymidineStates[5][i] = State( None, name='T_SE_pos_'+str(i) )                                                          #5 SE

 		#tie the emission probability of M1 to M2 so that they update together as the model learns
		thymidineStates[4][i].tie( thymidineStates[3][i] )                  

		#tie the insertion probabilities together so that we just get one uniform distribution
		if i != 0:
			thymidineStates[2][0].tie( thymidineStates[2][i] )

		#add the state objects to the hmm model object
		for j in range(6):
			hmm.add_state(thymidineStates[j][i])

		#internal transitions for the thymidine branch branch.  Group them so that they update together.
		#from SS
		hmm.add_transition(thymidineStates[0][i], thymidineStates[3][i], internalSS2M1 , group='internal_SS-to-M1')
		hmm.add_transition(thymidineStates[0][i], thymidineStates[4][i], internalSS2M2, group='internal_SS-to-M2')

		#from D
		hmm.add_transition(thymidineStates[1][i], thymidineStates[2][i], internalD2I, group='internal_D-to-I')

		#from I
		hmm.add_transition(thymidineStates[2][i], thymidineStates[2][i], internalI2I, group='internal_I-to-I')
		hmm.add_transition(thymidineStates[2][i], thymidineStates[0][i], internalI2SS, group='internal_I-to-SS')

		#from M1
		hmm.add_transition(thymidineStates[3][i], thymidineStates[3][i], internalM12M1, group='internal_M1-to-M1')
		hmm.add_transition(thymidineStates[3][i], thymidineStates[5][i], internalM12SE, group='internal_M1-to-SE')

		#from M2
		hmm.add_transition(thymidineStates[4][i], thymidineStates[4][i], internalM22M2, group='internal_M2-to-M2')
		hmm.add_transition(thymidineStates[4][i], thymidineStates[5][i], internalM22SE, group='internal_M2-to-SE')

		#from SE
		hmm.add_transition(thymidineStates[5][i], thymidineStates[2][i], internalSE2I, group='internal_SE-to-I')


		

	###########################################################################################################################
	# Handle Transitions Between Modules, Handle Analogue Branch

	#We have to reach forward to the next state (state i+1) so it's easier to just do this in a separate loop from the internal transitions one
	for i, char in enumerate(refSequence[:-5]): 

		###################################################
		# Handle Base Analogue Branch
		#    / [] - [] \
		# [] - [] - [] - []

		#If we're making an analogue branched HMM, there is a T about to enter position 3, 
		#we're more than 4 spaces from the end (we need to reach forward to i+3 to add the 
		#analogue branch), and we have data for at least one of the BrdU positions, then
		#we need a parallel analogue branch.

		#check if there is a T (which could be a base analogue) coming into position 4, and make sure we have analogue data for the analogue at position 3 and 4
		#if at least one of the modules in the branch could have a Hellinger distance of >0.5, make that branch
		makeFork = False
		if refSequence[i + 4] == 'T':
			analogue6merPos4 = refSequence[i+1:i+4] + 'B' + refSequence[i+5:i+7]
			analogue6merPos3 = refSequence[i+2:i+4] + 'B' + refSequence[i+5:i+8]
			if (analogue6merPos4 in poreModel) and (analogue6merPos3 in poreModel):
				makeFork = True
				analoguePositions.append(i + 1)
				analoguePositions.append(i + 2)

		#if we pass the check above and we're not at the end, create the forking structure for a base analogue 
		if makeFork and (i <= last_i-4):

			#assign a list of analogue state objects to the dictionary, keyed by the position of the 6mer in the reference
			analogueStates[i] = [	[State( None, name='B_SS_pos_'+str(i+1)+'_branchFrom'+str(i) ),                                     #0 SS
					     	State( None, name='B_D_pos_'+str(i+1)+'_branchFrom'+str(i) ),                                      #1 D
					     	State( UniformDistribution(30, 130, frozen=True), name='B_I_pos_'+str(i+1)+'_branchFrom'+str(i) ),               #2 I
					     	State( None, name='B_M1_pos_'+str(i+1)+'_branchFrom'+str(i) ),                                      #3 M1 - we'll tie this to M2 in a minute
					     	State( NormalDistribution( poreModel[analogue6merPos4][0], poreModel[analogue6merPos4][1] ), name='B_M2_pos_'+str(i+1)+'_branchFrom'+str(i) ),               #4 M2
					     	State( None, name='B_SE_pos_'+str(i+1)+'_branchFrom'+str(i) )],                                     #5 SE
						[State( None, name='B_SS_pos_'+str(i+2)+'_branchFrom'+str(i) ),                                     #0 SS
					     	State( None, name='B_D_pos_'+str(i+2)+'_branchFrom'+str(i) ),                                      #1 D
					     	State( UniformDistribution(30, 130, frozen=True), name='B_I_pos_'+str(i+2)+'_branchFrom'+str(i) ),               #2 I
					     	State( None, name='B_M1_pos_'+str(i+2)+'_branchFrom'+str(i) ),                                      #3 M1 - we'll tie this to M2 in a minute
					     	State( NormalDistribution( poreModel[analogue6merPos3][0], poreModel[analogue6merPos3][1] ), name='B_M2_pos_'+str(i+2)+'_branchFrom'+str(i) ),               #4 M2
					     	State( None, name='B_SE_pos_'+str(i+2)+'_branchFrom'+str(i) )]	]                             #5 SE
			
			#tie M1 to M2
			analogueStates[i][0][4].tie(analogueStates[i][0][3])
			analogueStates[i][1][4].tie(analogueStates[i][1][3])

			#tie the insertion probabilities together 
			thymidineStates[2][0].tie( analogueStates[i][0][2] )
			thymidineStates[2][0].tie( analogueStates[i][1][2] )

			#add the base analogue states to the hmm
			for j in range(6):
				hmm.add_states(analogueStates[i][0][j])
				hmm.add_states(analogueStates[i][1][j])

			
			#internal transitions both base analogue modules in the fork 
			for j in range(0,2):
				#from SS
				hmm.add_transition(analogueStates[i][j][0], analogueStates[i][j][3], internalSS2M1, group='internal_SS-to-M1')
				hmm.add_transition(analogueStates[i][j][0], analogueStates[i][j][4], internalSS2M2, group='internal_SS-to-M2')

				#from D
				hmm.add_transition(analogueStates[i][j][1], analogueStates[i][j][2], internalD2I, group='internal_D-to-I')

				#from I
				hmm.add_transition(analogueStates[i][j][2], analogueStates[i][j][2], internalI2I, group='internal_I-to-I')
				hmm.add_transition(analogueStates[i][j][2], analogueStates[i][j][0], internalI2SS, group='internal_I-to-SS')

				#from M1
				hmm.add_transition(analogueStates[i][j][3], analogueStates[i][j][3], internalM12M1, group='internal_M1-to-M1')
				hmm.add_transition(analogueStates[i][j][3], analogueStates[i][j][5], internalM12SE, group='internal_M1-to-SE')

				#from M2
				hmm.add_transition(analogueStates[i][j][4], analogueStates[i][j][4], internalM22M2, group='internal_M2-to-M2')
				hmm.add_transition(analogueStates[i][j][4], analogueStates[i][j][5], internalM22SE, group='internal_M2-to-SE')

				#from SE
				hmm.add_transition(analogueStates[i][j][5], analogueStates[i][j][2], internalSE2I, group='internal_SE-to-I')

			# create external module transitions for:
			# - leaving the thymidine branch at position i
			# - analogue-to-analogue in the analogue branch
			# - analogue at position i+2 merging back to the thymidine branch at thymidine module i+3

			#    */ []  - []  \
			# []  - []  - []  - []
			#T module to analogue module at position i+1
			hmm.add_transition(thymidineStates[1][i], analogueStates[i][0][1], analogueConc*externalD2D, group='external_T2B_D-to-D')
			hmm.add_transition(thymidineStates[1][i], analogueStates[i][0][0], analogueConc*externalD2SS, group='external_T2B_D-to-SS')
			hmm.add_transition(thymidineStates[2][i], analogueStates[i][0][0], analogueConc*externalI2SS, group='external_T2B_I-to-SS')
			hmm.add_transition(thymidineStates[5][i], analogueStates[i][0][0], analogueConc*externalSE2SS, group='external_T2B_SE-to-SS')
			hmm.add_transition(thymidineStates[5][i], analogueStates[i][0][1], analogueConc*externalSE2D, group='external_T2B_SE-to-D')

			#     / []  - []  \
			# [] *- []  - []  - []
			#T module to T module
			hmm.add_transition(thymidineStates[1][i], thymidineStates[1][i+1], (1-analogueConc)*externalD2D, group='external_T2T_D-to-D')
			hmm.add_transition(thymidineStates[1][i], thymidineStates[0][i+1], (1-analogueConc)*externalD2SS, group='external_T2T_D-to-SS')
			hmm.add_transition(thymidineStates[2][i], thymidineStates[0][i+1], (1-analogueConc)*externalI2SS, group='external_T2T_I-to-SS')
			hmm.add_transition(thymidineStates[5][i], thymidineStates[0][i+1], (1-analogueConc)*externalSE2SS, group='external_T2T_SE-to-SS')
			hmm.add_transition(thymidineStates[5][i], thymidineStates[1][i+1], (1-analogueConc)*externalSE2D, group='external_T2T_SE-to-D')

			#     / [] *- []  \
			# []  - []  - []  - []
			#analogue module at i+1 to analogue module at i+2
			hmm.add_transition(analogueStates[i][0][1], analogueStates[i][1][1], externalD2D, group='external_B2B_D-to-D')
			hmm.add_transition(analogueStates[i][0][1], analogueStates[i][1][0], externalD2SS, group='external_B2B_D-to-SS')
			hmm.add_transition(analogueStates[i][0][2], analogueStates[i][1][0], externalI2SS, group='external_B2B_I-to-SS')
			hmm.add_transition(analogueStates[i][0][5], analogueStates[i][1][0], externalSE2SS, group='external_B2B_SE-to-SS')
			hmm.add_transition(analogueStates[i][0][5], analogueStates[i][1][1], externalSE2D, group='external_B2B_SE-to-D')

			#     / []  - [] *\
			# []  - []  - []  - []
			hmm.add_transition(analogueStates[i][1][1], thymidineStates[1][i+3], externalD2D, group='external_B2T_D-to-D')
			hmm.add_transition(analogueStates[i][1][1], thymidineStates[0][i+3], externalD2SS, group='external_B2T_D-to-SS')
			hmm.add_transition(analogueStates[i][1][2], thymidineStates[0][i+3], externalI2SS, group='external_B2T_I-to-SS')
			hmm.add_transition(analogueStates[i][1][5], thymidineStates[0][i+3], externalSE2SS, group='external_B2T_SE-to-SS')
			hmm.add_transition(analogueStates[i][1][5], thymidineStates[1][i+3], externalSE2D, group='external_B2T_SE-to-D')

		#Otherwise (analogue = None or there is no T at position 3) just transition to the next T module
		#But don't execute this if we're at the end, because there's no i+1 to reach forward to.
		elif i != last_i:
			# [] *- []
			#T module to T module
			hmm.add_transition(thymidineStates[1][i], thymidineStates[1][i+1], externalD2D, group='external_T2T_D-to-D')
			hmm.add_transition(thymidineStates[1][i], thymidineStates[0][i+1], externalD2SS, group='external_T2T_D-to-SS')
			hmm.add_transition(thymidineStates[2][i], thymidineStates[0][i+1], externalI2SS, group='external_T2T_I-to-SS')
			hmm.add_transition(thymidineStates[5][i], thymidineStates[0][i+1], externalSE2SS, group='external_T2T_SE-to-SS')
			hmm.add_transition(thymidineStates[5][i], thymidineStates[1][i+1], externalSE2D, group='external_T2T_SE-to-D')


	###########################################################################################################################
	# Handle Start and End

	#add start states
	hmm.add_transition(hmm.start, thymidineStates[0][0], 0.5)
	hmm.add_transition(hmm.start, thymidineStates[1][0], 0.5)

	#handle end states
	hmm.add_transition(thymidineStates[1][last_i], hmm.end, externalD2D+externalD2SS)
	hmm.add_transition(thymidineStates[2][last_i], hmm.end, externalI2SS)
	hmm.add_transition(thymidineStates[5][last_i], hmm.end, externalSE2SS+externalSE2D)

	#bake the model
	t0 = time.time()
	hmm.bake(merge='None',verbose=True)
	t1 = time.time()
	print('Model optimised in '+str(t1-t0)+' seconds.')

	#uncomment if you want to autmoatically output the hmm object to a json file after it builds
	#with open('randIncHMM_untrained.json','w') as f:
	#	json.dump(hmm.to_json(),f)
	#f.close()

	return (analoguePositions, hmm)


def build_TrainingHMM(refSequence, poreModel):
#	Builds a HMM from a reference sequence with a topology appropriate for training analogue emission probability when the analogue is at a known, fixed position.  
#	Starts with the thymidine-containing kmer as a guess, and then refines the emission probability to that of the base analogue.
#	ARGUMENTS
#       ---------
#	- refSequence: reference sequence from import_reference
#	  type: string
#	- poreModel: the imported ONT pore model (from import_poreModel)
#	  type: dictionary
#	OUTPUTS
#       -------
#	- hmm: a training hidden Markov model for base analogues
#	  type: pomegranate object
	
	#set the hmm class
	hmm = HiddenMarkovModel() 

	refLength = len(refSequence)

	#two-dimensional array for the states, where the columns are the positions in the reference
	thymidineStates = [[0 for x in range(refLength)] for y in range(6)] 

	#Initial transitions within modules (internal transitions)
	internalSS2M1 = 0.97
	internalSS2M2 = 0.03
	internalD2I = 0.01
	internalI2I = 0.50
	internalI2SS = 0.49
	internalM12M1 = 0.51
	internalM12SE = 0.49
	internalM22M2 = 0.97
	internalM22SE = 0.03
	internalSE2I = 0.01

	#Initial transitions between modules (external transitions)
	externalD2D = 0.85
	externalD2SS = 0.14
	externalI2SS = 0.01
	externalSE2D = 0.12
	externalSE2SS = 0.87

	#this is really just a safety thing... get the last index in the iterator
	for last_i,x in enumerate(refSequence[:-5]): 
		pass
		
	###########################################################################################################################
	# Add States to Model
	
	#Create the HMM states.  Iterate through the reference sequence, and make the repeating HMM module for each position in the sequence.
	#Stop 7 characters before the end, because we'll need to reach forward by 7 to handle the branching.
	for i, char in enumerate(refSequence[:-5]): 

		###################################################
		# Handle Thymidine Branch States

		#grab the appropriate emission probabilities for the 6mer that we're on
		presentEmissions = poreModel[refSequence[i:i+6]]

		thymidineStates[0][i] = State( None, name='T_SS_pos_'+str(i) )                                                          #0 SS
		thymidineStates[1][i] = State( None, name='T_D_pos_'+str(i) )                                                           #1 D
		thymidineStates[2][i] = State( UniformDistribution(30, 130, frozen=True), name='T_I_pos_'+str(i) )                                   #2 I  - uniformly distributed on 30pA to 130pA
		thymidineStates[3][i] = State( None, name='T_M1_pos_'+str(i) )                                                          #3 M1 - we'll tie this to M2 in a minute
		thymidineStates[4][i] = State( NormalDistribution(presentEmissions[0], presentEmissions[1]), name='T_M2_pos_'+str(i) )  #4 M2 
		thymidineStates[5][i] = State( None, name='T_SE_pos_'+str(i) )                                                          #5 SE

 		#tie the emission probability of M1 to M2 so that they update together as the model learns
		thymidineStates[4][i].tie( thymidineStates[3][i] )                  

		#tie the insertion probabilities together so that we just get one uniform distribution
		if i != 0:
			thymidineStates[2][0].tie( thymidineStates[2][i] )

		#add the state objects to the hmm model object
		for j in range(6):
			hmm.add_state(thymidineStates[j][i])

		#internal transitions for the thymidine branch branch.  Group them so that they update together.
		#from SS
		hmm.add_transition(thymidineStates[0][i], thymidineStates[3][i], internalSS2M1 , group='internal_SS-to-M1')
		hmm.add_transition(thymidineStates[0][i], thymidineStates[4][i], internalSS2M2, group='internal_SS-to-M2')

		#from D
		hmm.add_transition(thymidineStates[1][i], thymidineStates[2][i], internalD2I, group='internal_D-to-I')

		#from I
		hmm.add_transition(thymidineStates[2][i], thymidineStates[2][i], internalI2I, group='internal_I-to-I')
		hmm.add_transition(thymidineStates[2][i], thymidineStates[0][i], internalI2SS, group='internal_I-to-SS')

		#from M1
		hmm.add_transition(thymidineStates[3][i], thymidineStates[3][i], internalM12M1, group='internal_M1-to-M1')
		hmm.add_transition(thymidineStates[3][i], thymidineStates[5][i], internalM12SE, group='internal_M1-to-SE')

		#from M2
		hmm.add_transition(thymidineStates[4][i], thymidineStates[4][i], internalM22M2, group='internal_M2-to-M2')
		hmm.add_transition(thymidineStates[4][i], thymidineStates[5][i], internalM22SE, group='internal_M2-to-SE')

		#from SE
		hmm.add_transition(thymidineStates[5][i], thymidineStates[2][i], internalSE2I, group='internal_SE-to-I')


		

	###########################################################################################################################
	# Handle Transitions Between Modules, Handle Analogue Branch

	#We have to reach forward to the next state (state i+1) so it's easier to just do this in a separate loop from the internal transitions one
	for i, char in enumerate(refSequence[:-5]): 

		#Don't execute this if we're at the end, because there's no i+1 to reach forward to.
		if i != last_i:
			# [] *- []
			#T module to T module
			hmm.add_transition(thymidineStates[1][i], thymidineStates[1][i+1], externalD2D, group='external_T2T_D-to-D')
			hmm.add_transition(thymidineStates[1][i], thymidineStates[0][i+1], externalD2SS, group='external_T2T_D-to-SS')
			hmm.add_transition(thymidineStates[2][i], thymidineStates[0][i+1], externalI2SS, group='external_T2T_I-to-SS')
			hmm.add_transition(thymidineStates[5][i], thymidineStates[0][i+1], externalSE2SS, group='external_T2T_SE-to-SS')
			hmm.add_transition(thymidineStates[5][i], thymidineStates[1][i+1], externalSE2D, group='external_T2T_SE-to-D')


	###########################################################################################################################
	# Handle Start and End

	#add start states
	hmm.add_transition(hmm.start, thymidineStates[0][0], 0.5)
	hmm.add_transition(hmm.start, thymidineStates[1][0], 0.5)

	#handle end states
	hmm.add_transition(thymidineStates[1][last_i], hmm.end, externalD2D+externalD2SS)
	hmm.add_transition(thymidineStates[2][last_i], hmm.end, externalI2SS)
	hmm.add_transition(thymidineStates[5][last_i], hmm.end, externalSE2SS+externalSE2D)

	#bake the model
	t0 = time.time()
	hmm.bake(merge='None',verbose=True)
	t1 = time.time()
	print('Model optimised in '+str(t1-t0)+' seconds.')

	return hmm
	

def build_softHMM(refSequence, poreModel, analogueEmissions=None, analoguePos=None):
#	Builds a HMM from a reference sequence with a topology appropriate for training analogue emission probability when the analogue is at a known, fixed position.  
#	Starts with the thymidine-containing kmer as a guess, and then refines the emission probability to that of the base analogue.
#	ARGUMENTS
#       ---------
#	- refSequence: reference sequence from import_reference
#	  type: string
#	- poreModel: the imported ONT pore model (from import_poreModel)
#	  type: dictionary
#	OUTPUTS
#       -------
#	- hmm: a training hidden Markov model for base analogues
#	  type: pomegranate object
	
	#set the hmm class
	hmm = HiddenMarkovModel() 

	refLength = len(refSequence)

	#two-dimensional array for the states, where the columns are the positions in the reference
	thymidineStates = [[0 for x in range(refLength)] for y in range(6)] 

	#Initial transitions within modules (internal transitions)
	internalSS2M1 = 0.97
	internalSS2M2 = 0.03
	internalD2I = 0.01
	internalI2I = 0.50
	internalI2SS = 0.49
	internalM12M1 = 0.51
	internalM12SE = 0.49
	internalM22M2 = 0.97
	internalM22SE = 0.03
	internalSE2I = 0.01

	#Initial transitions between modules (external transitions)
	externalD2D = 0.85
	externalD2SS = 0.14
	externalI2SS = 0.01
	externalSE2D = 0.12
	externalSE2SS = 0.87

	#this is really just a safety thing... get the last index in the iterator
	for last_i,x in enumerate(refSequence[:-5]): 
		pass
		
	###########################################################################################################################
	# Add States to Model
	
	#Create the HMM states.  Iterate through the reference sequence, and make the repeating HMM module for each position in the sequence.
	#Stop 7 characters before the end, because we'll need to reach forward by 7 to handle the branching.
	for i, char in enumerate(refSequence[:-5]): 

		###################################################
		# Handle Thymidine Branch States

		#grab the appropriate emission probabilities for the 6mer that we're on
		if analogueEmissions is not None:
			if i == analougePos - 3:
				presentEmissions = analogueEmissions[ reference[ i : analoguePos ] + 'B' + reference[ analogue + 1 : analogue + 3 ] ]
			elif i == analoguePos - 2:
				presentEmissions = analogueEmissions[ reference[ i : analoguePos ] + 'B' + reference[ analogue + 1 : analogue + 4 ] ]				 
		else:
			presentEmissions = poreModel[refSequence[i:i+6]]

		thymidineStates[0][i] = State( None, name='T_SS_pos_'+str(i) )                                                          #0 SS
		thymidineStates[1][i] = State( None, name='T_D_pos_'+str(i) )                                                           #1 D
		thymidineStates[2][i] = State( UniformDistribution(30, 130, frozen=True), name='T_I_pos_'+str(i) )                                   #2 I  - uniformly distributed on 30pA to 130pA
		thymidineStates[3][i] = State( None, name='T_M1_pos_'+str(i) )                                                          #3 M1 - we'll tie this to M2 in a minute
		thymidineStates[4][i] = State( NormalDistribution(presentEmissions[0], presentEmissions[1]), name='T_M2_pos_'+str(i) )  #4 M2 
		thymidineStates[5][i] = State( None, name='T_SE_pos_'+str(i) )                                                          #5 SE

 		#tie the emission probability of M1 to M2 so that they update together as the model learns
		thymidineStates[4][i].tie( thymidineStates[3][i] )                  

		#tie the insertion probabilities together so that we just get one uniform distribution
		if i != 0:
			thymidineStates[2][0].tie( thymidineStates[2][i] )

		#add the state objects to the hmm model object
		for j in range(6):
			hmm.add_state(thymidineStates[j][i])

		#internal transitions for the thymidine branch branch.  Group them so that they update together.
		#from SS
		hmm.add_transition(thymidineStates[0][i], thymidineStates[3][i], internalSS2M1 , group='internal_SS-to-M1')
		hmm.add_transition(thymidineStates[0][i], thymidineStates[4][i], internalSS2M2, group='internal_SS-to-M2')

		#from D
		hmm.add_transition(thymidineStates[1][i], thymidineStates[2][i], internalD2I, group='internal_D-to-I')

		#from I
		hmm.add_transition(thymidineStates[2][i], thymidineStates[2][i], internalI2I, group='internal_I-to-I')
		hmm.add_transition(thymidineStates[2][i], thymidineStates[0][i], internalI2SS, group='internal_I-to-SS')

		#from M1
		hmm.add_transition(thymidineStates[3][i], thymidineStates[3][i], internalM12M1, group='internal_M1-to-M1')
		hmm.add_transition(thymidineStates[3][i], thymidineStates[5][i], internalM12SE, group='internal_M1-to-SE')

		#from M2
		hmm.add_transition(thymidineStates[4][i], thymidineStates[4][i], internalM22M2, group='internal_M2-to-M2')
		hmm.add_transition(thymidineStates[4][i], thymidineStates[5][i], internalM22SE, group='internal_M2-to-SE')

		#from SE
		hmm.add_transition(thymidineStates[5][i], thymidineStates[2][i], internalSE2I, group='internal_SE-to-I')

	###########################################################################################################################
	# Handle Transitions Between Modules, Handle Analogue Branch

	#We have to reach forward to the next state (state i+1) so it's easier to just do this in a separate loop from the internal transitions one
	for i, char in enumerate(refSequence[:-5]): 

		#Don't execute this if we're at the end, because there's no i+1 to reach forward to.
		if i != last_i:
			# [] *- []
			#T module to T module
			hmm.add_transition(thymidineStates[1][i], thymidineStates[1][i+1], externalD2D, group='external_T2T_D-to-D')
			hmm.add_transition(thymidineStates[1][i], thymidineStates[0][i+1], externalD2SS, group='external_T2T_D-to-SS')
			hmm.add_transition(thymidineStates[2][i], thymidineStates[0][i+1], externalI2SS, group='external_T2T_I-to-SS')
			hmm.add_transition(thymidineStates[5][i], thymidineStates[0][i+1], externalSE2SS, group='external_T2T_SE-to-SS')
			hmm.add_transition(thymidineStates[5][i], thymidineStates[1][i+1], externalSE2D, group='external_T2T_SE-to-D')


	###########################################################################################################################
	# Handle Start and End

	#add beginning and end garbage collection states
	gcStart = State( UniformDistribution(30, 130, frozen=True), name='startGarbageCollection' )
	gcEnd = State( UniformDistribution(30, 130, frozen=True), name='endGarbageCollection' )
	hmm.add_state( gcStart )
	hmm.add_state( gcEnd )

	#garbage collection start transitions
	hmm.add_transition(hmm.start, gcStart, 1.0 )
	hmm.add_transition(gcStart, gcStart, 0.8 )
	hmm.add_transition(gcStart, thymidineStates[0][0], 0.1)
	hmm.add_transition(gcStart, thymidineStates[1][0], 0.1)

	#handle garbage collection end states
	hmm.add_transition(gcEnd, gcEnd, 0.8 )
	hmm.add_transition(gcEnd, hmm.end, 0.2 )
	hmm.add_transition(thymidineStates[1][last_i], gcEnd, externalD2D+externalD2SS)
	hmm.add_transition(thymidineStates[2][last_i], gcEnd, externalI2SS)
	hmm.add_transition(thymidineStates[5][last_i], gcEnd, externalSE2SS+externalSE2D)

	#bake the model
	t0 = time.time()
	hmm.bake(merge='None',verbose=True)
	t1 = time.time()
	print('Model optimised in '+str(t1-t0)+' seconds.')

	return hmm

		
		
