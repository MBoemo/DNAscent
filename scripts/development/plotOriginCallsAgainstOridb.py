import sys
import matplotlib
matplotlib.use('Agg')
from matplotlib import pyplot as plt
import numpy as np
import random
import scipy
from scipy import special



def foldedMean(oriMean, oriStd):

	return oriStd*np.sqrt( 2.0 / np.pi ) * np.exp( - oriMean**2 / (2*oriStd**2) ) + oriMean * special.erf( oriMean / np.sqrt( 2*oriStd**2 ) )

def foldedStd(oriMean, oriStd, fm):

	return np.sqrt( oriMean**2 + oriStd**2 - fm**2 )

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
f = open(sys.argv[1],'r')
calledOriPos = {}
calledOriPosWithDistance = {}
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
		calledOriPosWithDistance[chromosome] = []
f.close()

#get the origins from Oridb
f = open(sys.argv[2],'r')
dbOriPos = {}
for line in f:

	splitLine = line.rstrip().split('\t')
	chromosome = 'chr'+int2RomanNumeral(int(splitLine[0][3:]))
	loc = (int(splitLine[1]), int(splitLine[2]))

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

			distFromOri = min( abs(oriCall - knownOri[0]), abs(oriCall - knownOri[1]) ) #0 if within bounds

			if distFromOri < abs(minDist):

				if oriCall > knownOri[0] and oriCall < knownOri[1]:

					minDist = 0

				else:
					if abs(oriCall - knownOri[0]) < abs(oriCall - knownOri[1]):

						minDist = oriCall - knownOri[0]
					else:
						minDist = oriCall - knownOri[1]

		distances.append(minDist)
		calledOriPosWithDistance[chromosome].append( (oriCall, minDist) )


	#randomly generated
	for i in range(0, len(calledOriPos[chromosome])):
		randomPos = random.randint(0,getChromosomeLength(chromosome))
		minDist = 100000000
		for knownOri in dbOriPos[chromosome]:

			distFromOri = min( abs(randomPos - knownOri[0]), abs(randomPos - knownOri[1]) )

			if distFromOri < abs(minDist):

				if randomPos > knownOri[0] and randomPos < knownOri[1]:

					minDist = 0

				else:
					if abs(randomPos - knownOri[0]) < abs(randomPos - knownOri[1]):

						minDist = randomPos - knownOri[0]
					else:
						minDist = randomPos - knownOri[1]

		randomDistances.append(minDist)

	
	

plt.figure(1)
plt.hist(distances, 150,color='k',label="Called by Osiris",linewidth=0,alpha=0.3)
plt.hist(randomDistances, 150,color='r',label="Random",linewidth=0,alpha=0.3)
plt.legend()
plt.xlabel('Distance to Closest Origin (bp)')
plt.ylabel('Count')
plt.xlim(-20000,20000)
plt.savefig('originCallDistance.pdf')
plt.cla()

plt.figure(2)
plt.hist(np.absolute(distances), 50,color='k',label="Called by Osiris",linewidth=0,alpha=0.3)
plt.hist(np.absolute(randomDistances), 50,color='r',label="Random",linewidth=0,alpha=0.3)
plt.legend()
plt.xlabel('Distance to Closest Origin (bp)')
plt.ylabel('Count')
plt.xlim(0,25000)
plt.savefig('originCallDistance_folded.pdf')
plt.cla()


oriMean = np.mean(distances)
oriStd = np.std(distances)

#compute folded normal
fm = oriStd*np.sqrt( 2.0 / np.pi ) * np.exp( - oriMean**2 / (2*oriStd**2) ) + oriMean * special.erf( oriMean / np.sqrt( 2*oriStd**2 ) )
fstd = np.sqrt( oriMean**2 + oriStd**2 - fm**2 )
print "Osiris folded mean, stdv:",fm,fstd
print "Osiris normal mean, stdv:",oriMean, oriStd
print "Osiris abs mean, stdv:",np.mean(np.absolute(distances)),np.std(np.absolute(distances))
print "Osiris median:",np.median(np.absolute(distances))
print " "
randomMean = np.mean(randomDistances)
randomStd = np.std(randomDistances)
randomFoldedMean = foldedMean(randomMean, randomStd)
print "Random normal mean, stdv:", randomMean, randomStd
print "Random folded mean, stdv:", randomFoldedMean, foldedStd(randomMean,randomStd,randomFoldedMean)
print "Random median:",np.median(np.absolute(randomDistances))

maxWindowSize = 11000
closeMinDistances = []
farMinDistances = []
ynIsolatedFar = [0.0] * len(range(1000,maxWindowSize,1000))
ynIsolatedClose = [0.0] * len(range(1000,maxWindowSize,1000))
ydIsolatedFar = [0.0] * len(range(1000,maxWindowSize,1000))
ydIsolatedClose = [0.0] * len(range(1000,maxWindowSize,1000))

for chromosome in calledOriPosWithDistance:

	close = []
	far = []

	for oriDistPair in calledOriPosWithDistance[chromosome]:

		if (abs(oriDistPair[1]) > fm):
		#if abs(oriDistPair[1]) > np.median(np.absolute(distances)):

			far.append(oriDistPair[0])
		else:
			close.append(oriDistPair[0])

	for idx, windowSize in enumerate(range(1000,maxWindowSize,1000)):

		isolatedClose = 0
		isolatedFar = 0

		for i, c in enumerate(close):
			closeMinDistance = 10000000000
			partnerFound = False
			for d in close:
			
				if d == c:
					continue

				if abs(c-d) < windowSize:
					closeMinDistance = abs(c-d)
					partnerFound = True
			for d in far:

				if abs(c-d) < windowSize:
					closeMinDistance = abs(c-d)
					partnerFound = True


			if partnerFound:
				closeMinDistances.append( closeMinDistance )
			else:
				isolatedClose += 1


		for i, c in enumerate(far):
			farMinDistance = 10000000000
			partnerFound = False
			for d in far:

				if d == c:
					continue

				if abs(c-d) < windowSize:
					farMinDistance = abs(c-d)
					partnerFound = True
			for d in close:

				if abs(c-d) < windowSize:
					farMinDistance = abs(c-d)
					partnerFound = True

			if partnerFound:
				farMinDistances.append( farMinDistance )
			else: 
				isolatedFar += 1

		ynIsolatedFar[idx] += float(isolatedFar)
		ydIsolatedFar[idx] += float(len(far))
		ynIsolatedClose[idx] += float(isolatedClose)
		ydIsolatedClose[idx] +=float(len(close))

print "number of close:",ydIsolatedClose[0]
print "number of far:",ydIsolatedFar[0]
plt.figure(3)
plt.plot( range(1000,maxWindowSize,1000), np.divide(np.array(ynIsolatedClose), np.array(ydIsolatedClose)), label = 'Close to Known Origins')
plt.plot( range(1000,maxWindowSize,1000), np.divide(np.array(ynIsolatedFar), np.array(ydIsolatedFar)), label = 'Far from Known Origins')
plt.legend()
plt.xlabel('Window Size (bp)')
plt.ylabel('Isolated Fraction')
plt.savefig('originCallDistance_isolatedByWindow.pdf')


#fig, ax1 = plt.subplots()
#ax1.hist(closeMinDistances,100,color='b',linewidth=0,alpha=0.3)
#ax1.set_xlabel('Distance From Closest Origin (bp)')
#ax1.set_ylabel('Inside One Stdv', color='b')
#ax1.tick_params('y', colors='b')

#ax2 = ax1.twinx()
#ax2.hist(farMinDistances,100,color='r',linewidth=0,alpha=0.3)
#ax2.set_ylabel('Outside One Stdv', color='r')
#ax2.tick_params('y', colors='r')
#plt.xlim(0,4000)
#plt.savefig('originCallDistance_nearAndFar.pdf')


	

		
