import sys
import itertools
import random
import numpy as np
np.set_printoptions(threshold=sys.maxsize)
import sys
import os
import pickle
import random
import math

barcode = int(sys.argv[1])

bc08folderPath = '/home/mb915/rds/rds-mb915-notbackedup/data/2018_06_18_CAM_ONT_gDNA_BrdU_40_60_80_100_full/cnn_training/data_8features_bc08_withStd'
bc10folderPath = '/home/mb915/rds/rds-mb915-notbackedup/data/2018_06_18_CAM_ONT_gDNA_BrdU_40_60_80_100_full/cnn_training/data_8features_bc10_withStd'
bc11folderPath = '/home/mb915/rds/rds-mb915-notbackedup/data/2018_06_18_CAM_ONT_gDNA_BrdU_40_60_80_100_full/cnn_training/data_8features_bc11_withStd'
bc12folderPath = '/home/mb915/rds/rds-mb915-notbackedup/data/2018_06_18_CAM_ONT_gDNA_BrdU_40_60_80_100_full/cnn_training/data_8features_bc12_withStd'

bc08dnascent = '/home/mb915/rds/rds-mb915-notbackedup/data/2018_06_18_CAM_ONT_gDNA_BrdU_40_60_80_100_full/barcode08/commit620d798_l_4000_q_20.barcode08.trainingData'
bc10dnascent = '/home/mb915/rds/rds-mb915-notbackedup/data/2018_06_18_CAM_ONT_gDNA_BrdU_40_60_80_100_full/barcode10/commit620d798_l_4000_q_20.trainingData'
bc11dnascent = '/home/mb915/rds/rds-mb915-notbackedup/data/2018_06_18_CAM_ONT_gDNA_BrdU_40_60_80_100_full/barcode11/commit620d798_l_4000_q_20.barcode11.trainingData'
bc12dnascent = '/home/mb915/rds/rds-mb915-notbackedup/data/2018_06_18_CAM_ONT_gDNA_BrdU_40_60_80_100_full/barcode12/commit620d798_l_4000_q_20.trainingData'

inputFiles = [(0., bc08dnascent, bc08folderPath),(0.26, bc10dnascent, bc10folderPath),(0.5, bc11dnascent, bc11folderPath),(0.8, bc12dnascent, bc12folderPath)]
maxLen = 4000
maxReads = int(sys.argv[2])

llThreshold = 1.25


#-------------------------------------------------
#
class trainingRead:

	def __init__(self, read_sixMers, read_eventMeans, read_eventStd, read_eventLength, read_modelMeans, read_modelStd, read_positions, readID, analogueConc, logLikelihood):
		self.sixMers = read_sixMers[0:maxLen]
		self.eventMean = read_eventMeans[0:maxLen]
		self.eventStd = read_eventStd[0:maxLen]
		self.eventLength = read_eventLength[0:maxLen]
		self.modelMeans = read_modelMeans[0:maxLen]
		self.modelStd = read_modelStd[0:maxLen]
		self.logLikelihood = logLikelihood[0:maxLen]
		self.readID = readID
		self.analogueConc = analogueConc

		gaps = [1]
		for i in range(1, maxLen):
			gaps.append( abs(read_positions[i] - read_positions[i-1]) )
		self.gaps = gaps

		'''
		print(self.sixMers[0:5])
		print(self.eventMean[0:5])
		print(self.eventStd[0:5])
		print(self.stutter[0:5])
		print(self.eventLength[0:5])
		print(self.modelMeans[0:5])
		print(self.modelStd[0:5])
		print(self.gaps[0:5])
		print(read_positions[0:5])
		print(logLikelihood[0:5])
		print('-----------------------')
		'''

		if not len(self.sixMers) == len(self.eventMean) == len(self.eventStd) == len(self.eventLength) == len(self.modelMeans) == len(self.modelStd) == len(self.gaps) == len(self.logLikelihood):
			print(len(self.sixMers), len(self.logLikelihood))
			print("Length Mismatch")
			sys.exit()


#-------------------------------------------------
#one-hot for bases
baseToInt = {'A':0, 'T':1, 'G':2, 'C':3}
intToBase = {0:'A', 1:'T', 2:'G', 3:'C'}


#-------------------------------------------------
#reshape into tensors
def saveRead(trainingRead, readID, folderPath):

	f = open(folderPath + '/' + readID + '.p', 'wb')
	pickle.dump(trainingRead, f, protocol=pickle.HIGHEST_PROTOCOL)
	f.close()


#-------------------------------------------------
#reshape into tensors
def weighted_avg_and_std(values, weights):
	"""
	Return the weighted average and standard deviation.

	values, weights -- Numpy ndarrays with the same shape.
	"""

	values = np.array(values)
	weights = np.array(weights)

	average = np.average(values, weights=weights)
	# Fast and numerically precise:
	variance = np.average((values-average)**2, weights=weights)
	return math.sqrt(variance)


#-------------------------------------------------
#pull from training data file
def importFromFile(analogueConc, fname, folderPath):

	#count the number of reads in order to randomize across the genome
	totalReads = 0
	readsToUse = []	
	f = open(fname,'r')
	for line in f:
		if line[0] == '>':
			totalReads += 1
	f.close()
	readsToUse = random.sample(range(0,totalReads), maxReads)


	print('Parsing training data file...')

	read_sixMers = []
	read_eventMeans = []
	read_eventStd = []
	read_stutter = []
	read_lengthMeans = []
	read_lengthStd = []
	read_modelMeans = []
	read_modelStd = []
	read_positions = []
	read_BrdUcalls = []

	pos_eventMeans = []
	pos_lengths = []

	prevReadID = ''
	prevPos = -1
	
	readsLoaded = 0
	switch = False
	f = open(fname,'r')
	for line in f:

		if line[0] == '#':
			continue

		if line[0] == '>':

			readsLoaded += 1
			#print('Reads loaded: ',readsLoaded)

			#make an object out of this read
			if len(read_sixMers) > maxLen:

				if switch:
					tr = trainingRead(read_sixMers, read_eventMeans, read_eventStd, read_eventLength, read_modelMeans, read_modelStd, read_positions, prevReadID, analogueConc, read_BrdUcalls)
					saveRead(tr, prevReadID, folderPath)

			#reset for read
			read_sixMers = []
			read_eventMeans = []
			read_eventStd = []
			read_eventLength = []
			read_modelMeans = []
			read_modelStd = []
			read_positions = []
			read_BrdUcalls = []

			#reset for position
			pos_eventMeans = []
			pos_lengths = []
			prevPos = -1

			splitLine = line.rstrip().split()
			readID = splitLine[0][1:]
			chromosome = splitLine[1]
			mappingStart = int(splitLine[2])
			mappingEnd = int(splitLine[3])
			prevPos = -1

			if readsLoaded in readsToUse and chromosome != 'chrM':
				switch = True
			else:
				switch = False

		elif switch:

			splitLine = line.rstrip().split('\t')

			#skip insertion events
			if splitLine[4] == "NNNNNN":
				continue

			pos = int(splitLine[0])

			#don't get too close to the ends of the read, because we can't get an HMM detect there
			if abs(pos-mappingStart) < 26 or abs(pos-mappingEnd) < 26:
				continue

			#we've changed position
			if pos != prevPos and prevPos != -1:

				read_sixMers.append(sixMer)
				read_eventMeans.append(np.average(pos_eventMeans, weights=pos_lengths))
				read_eventStd.append(weighted_avg_and_std(pos_eventMeans, pos_lengths))
				read_eventLength.append(sum(pos_lengths))
				read_modelMeans.append(modelMean)
				read_modelStd.append(modelStd)
				read_positions.append(prevPos)
				read_BrdUcalls.append(BrdUcall)

				#reset for position
				pos_eventMeans = []
				pos_lengths = []
				prevPos = pos

			sixMer = splitLine[4]
			modelMean = float(splitLine[5])
			modelStd = float(splitLine[6])
			eventMean = float(splitLine[2])
			eventLength = float(splitLine[3])
			prevReadID = readID
			prevPos = pos

			#sort out BrdU calls
			BrdUcall = '-'
			if len(splitLine) == 8:
				BrdUcall = splitLine[7]
				
			pos_eventMeans.append(eventMean)
			pos_lengths.append(eventLength)

	f.close()
	print('Done.')


#-------------------------------------------------
#MAIN
tup = inputFiles[barcode]
importFromFile(tup[0], tup[1], tup[2])