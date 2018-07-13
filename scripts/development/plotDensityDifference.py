import sys
import matplotlib
matplotlib.use('Agg')
from matplotlib import pyplot as plt
import numpy as np
import matplotlib.mlab as mlab
from scipy import stats

#first argument is the workingData.osiris and second argument is the 5mer pore model

reference = 'AAGGTTAACCTGGTAACTGGGACACAAGACTCCAGCACCTGTAAAACGACGGCCAGTGAATTGTAATACGACTCACTATAGGGCGAATTGGGCCCGACGTCGCATGCTCCCGGCCGCCATGGCGGCCGCGGGAATTCGATTAAGTCAGACTCCCTTACAACACTTATGGTGAGACGACAGTTGGGTGGCCGCGCCTCATTGGAGATAACAGCCATCCTGATCCGCTGCGGTAAAACGTATATTCTCGATCTTTAACCAGTAGGTTTGGCAGTGAAGTTAGCAAGTACCCCTATGAACAATTCAAATGGGGACAAAATCCATGCTCTGTACGGAAGAGTTCCTAGCGCAAAGGAGGACGGCACTACATTATAGCTGGAATGCCTAAGCGACGCGAACCGAGGGTCTTGATACACGTCACACATGATGACATATCCCAGATTCGGGAAATAGTTTATTGAGTGGACCTGGCGAGCCGGGCGGGGGCTACCTTCGTAGATGTTTCTTAATCGTGCGTGGTATTACGCTCAGTCCGATAGACACCGGAGCTTTCGACCGTTGACCAAGCCTTGTGGGCAATCACGGGTTGCATCGCATACTAATTTAGAGAGGTGCTTCTAGTCGGCGTTACTCGTTCTGCACAGTATCACTAGTGAATTCGCGGCCGCCTGCAGGTCGACCATATGGGAGAGCTCCCAACGCGTTGGATGCATAGCTTGAGTATTCTATAGTGTCACCTAAATAGCTTGGCGTAATCATGGTCATAGCTGTTTCCTG'

#BrdU event alignment
f = open(sys.argv[1],'r')
pos2eventsBrdU = {}
for line in f:

	splitLine = line.rstrip().split(' ')
	position = splitLine[0]
	if position not in pos2eventsBrdU:
		pos2eventsBrdU[position] = []
	else:
		pos2eventsBrdU[position] += map(float,splitLine[1:])
f.close()


#thymidine event alignment
f = open(sys.argv[2],'r')
pos2eventsT = {}
for line in f:

	splitLine = line.rstrip().split(' ')
	position = splitLine[0]
	if position not in pos2eventsT:
		pos2eventsT[position] = []
	else:
		pos2eventsT[position] += map(float,splitLine[1:])
f.close()

#ONT model
model = {}
f = open(sys.argv[3],'r')
for line in f:

	if line[0] == '#' or line[0] == 'k':
		continue
	
	splitLine = line.rstrip().split('\t')
	model[splitLine[0]] = [	float(splitLine[1]), float(splitLine[2]) ]
f.close()

#plots
for i, key in enumerate(pos2eventsBrdU):
	
	fiveMer = reference[int(key):int(key)+5]
	x = np.linspace( np.mean(pos2eventsBrdU[key])-15, np.mean(pos2eventsBrdU[key])+15, 1000 )

	plt.figure(i)
	densityBrdU = stats.kde.gaussian_kde( pos2eventsBrdU[key] )
	densityT = stats.kde.gaussian_kde( pos2eventsT[key] )
	plt.plot( x, densityBrdU(x), label='BrdU Density')
	plt.plot( x, densityT(x), label='Thymidine Density')
	yModel = mlab.normpdf( x, model[fiveMer][0], model[fiveMer][1] )
	plt.plot( x, yModel, label='5mer Pore Model')

	densityDiff = np.array(densityBrdU(x) - densityT(x))
	densityDiff[ densityDiff < 0 ] = 0.0
	plt.plot( x, densityDiff, label='Density Difference')

	plt.xlabel('pA')
	plt.ylabel('Count')
	plt.title( reference[int(key):int(key)+5] )
	plt.legend(loc='upper right')
	plt.savefig( fiveMer + '.png' )
	plt.close()
