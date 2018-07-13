import sys
import matplotlib
matplotlib.use('Agg')
from matplotlib import pyplot as plt
import numpy as np
import matplotlib.mlab as mlab
from scipy import stats



#first argument is the workingData.osiris and second argument is the 5mer pore model


reference = 'AAGGTTAAAAGGTTACACAAACCCTGGACAAGCAGCACCTGTAAAACGACGGCCAGTGAATTGTAATACGACTCACTATAGGGCGAATTGGGCCCGACGTCGCATGCTCCCGGCCGCCATGGCGGCCGCGGGAATTCGATTATACTGTATTGATCGAGTTTGTGTTCTAATTGCGATACCATGATGGCTGAATAACCACCGAGGCAGACTAAGCCCGATGCAACGCATGTTTCGCGGATGAGTCCGTTAGGGGTGTCCTATAAGATATGTCACACTCCGGGACGAAGGTCGGCACCTCACGGGGCGGGTCTCAGGCGCGTACAACAGGAGCGCAGGTTCCCTGGTCAGTCAAGACGCCGGTTTTAAGGCTAGGTAGTGCGGCCTACTTACTATCCTCCCAAGGAATCGTTCATAGACAATCAGAATTTGAGCATTGGATTTCTTCCGAACTTGTTACGGCTCGCCAGTTGAAAGTGATAATGTGGCAAGCGGTCCATAAATATACGTGTAGATTACCGTCGCTGTGCTTTTGCTATTAACTAGAGTACCCAGATGTAAAGAGGGTTGGTGAAGCTACGAGAGAAGTCGTGGGGGATCTCCACTGCACGCTAAACGTCATTACTTTTTCAGCAGGCCCTGTGTATCACTAGTGAATTCGCGGCCGCCTGCAGGTCGACCATATGGGAGAGCTCCCAACGCGTTGGATGCATAGCTTGAGTATTCTATAGTGTCACCTAAATAGCTTGGCGTAATCATGGTCATAGCTGTTTCCTG'

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


model = {}
f = open(sys.argv[3],'r')
for line in f:

	if line[0] == '#' or line[0] == 'k':
		continue
	
	splitLine = line.rstrip().split('\t')
	model[splitLine[0]] = [	float(splitLine[1]), float(splitLine[2]) ]
f.close()

for i, key in enumerate(pos2eventsBrdU):
	
	sixMer = reference[int(key):int(key)+6]
	x = np.linspace( np.mean(pos2eventsBrdU[key])-15, np.mean(pos2eventsBrdU[key])+15, 1000 )

	plt.figure(i)
	densityBrdU = stats.kde.gaussian_kde( pos2eventsBrdU[key] )
	densityT = stats.kde.gaussian_kde( pos2eventsT[key] )
	plt.plot( x, densityBrdU(x), label='BrdU Density')
	plt.plot( x, densityT(x), label='Thymidine Density')
	yModel = mlab.normpdf( x, model[sixMer][0], model[sixMer][1] )
	plt.plot( x, yModel, label='6mer Pore Model')
	plt.xlabel('pA')
	plt.ylabel('Count')
	plt.title( reference[int(key):int(key)+6] )
	plt.legend(loc='upper right')
	plt.savefig( sixMer + '.png' )
	plt.close()
