import sys
import matplotlib
matplotlib.use('Agg')
from matplotlib import pyplot as plt
import numpy as np
import random

f = open(sys.argv[1],'r')

def int2RomanNumeral( number ):
	if number == 1:
		return "I"
	elif number == 2:
		return "II"
	elif number == 3:
		return "III"
	elif number == 4:
		return "IV"
	elif number == 5:
		return "V"
	elif number == 6:
		return "VI"
	elif number == 7:
		return "VII"
	elif number == 8:
		return "VIII"
	elif number == 9:
		return "IX"
	elif number == 10:
		return "X"
	elif number == 11:
		return "XI"
	elif number == 12:
		return "XII"
	elif number == 13:
		return "XIII"
	elif number == 14:
		return "XIV"
	elif number == 15:
		return "XV"
	elif number == 16:
		return "XVI"

def getChromosomeLength( chromosome ):
	
	if chromosome == "chrI":
		return 230218
	elif chromosome == "chrII":
		return 813184
	elif chromosome == "chrIII":
		return 316620
	elif chromosome == "chrIV":
		return 1531933
	elif chromosome == "chrV":
		return 576874
	elif chromosome == "chrVI":
		return 270161
	elif chromosome == "chrVII":
		return 1090940
	elif chromosome == "chrVIII":
		return 562643
	elif chromosome == "chrIX":
		return 439888
	elif chromosome == "chrX":
		return 745751
	elif chromosome == "chrXI":
		return 666816
	elif chromosome == "chrXII":
		return 1078177
	elif chromosome == "chrXIII":
		return 924431
	elif chromosome == "chrXIV":
		return 784333
	elif chromosome == "chrXV":
		return 1091291
	elif chromosome == "chrXVI":
		return 948066

#get the origin calls from Osiris
calledOriPos = {}
for line in f:
	if line[0] == '>':
		splitLine = line.rstrip().split(' ')
		splitLine2 = splitLine[1].split(':')
		chromosome = splitLine2[0]
	else:
		oriPos = int(line.rstrip())

		if chromosome in calledOriPos:
			calledOriPos[chromosome].append(oriPos)
		else:
			calledOriPos[chromosome] = [oriPos]
f.close()

#get the origins from Oridb
f = open(sys.argv[2],'r')
dbOriPos = {}
for line in f:
	splitLine = line.rstrip().split('\t')
	chromosome = 'chr'+int2RomanNumeral(int(splitLine[0][3:]))
	loc = (int(splitLine[1]) + int(splitLine[2]))/2
	if chromosome in dbOriPos:
		dbOriPos[chromosome].append(loc)
	else:
		dbOriPos[chromosome] = [loc]

#calculate nearest distances
distances = []
randomDistances = []
for chromosome in calledOriPos:

	#from Osiris
	for oriCall in calledOriPos[chromosome]:

		minDist = 100000000
		for knownOri in dbOriPos[chromosome]:
			if abs(oriCall - knownOri) < abs(minDist):
				minDist = oriCall - knownOri
		distances.append(minDist)

	#randomly generated
	for i in range(0, len(calledOriPos[chromosome])):
		randomPos = random.randint(0,getChromosomeLength(chromosome))
		minDist = 100000000
		for knownOri in dbOriPos[chromosome]:
			if abs(randomPos - knownOri) < abs(minDist):#
				minDist = randomPos - knownOri#
		randomDistances.append(minDist)

	
	

plt.figure()
plt.hist(distances, 100,color='k',label="Called by Osiris",linewidth=0,alpha=0.3)
plt.hist(randomDistances, 100,color='r',label="Random",linewidth=0,alpha=0.3)
plt.legend()
plt.xlabel('Distance to Closest Origin (bp)')
plt.ylabel('Count')
plt.savefig('originCallDistance.pdf')
print np.mean(distances)
print np.std(distances)
print np.mean(randomDistances)
print np.std(randomDistances)

	

		
