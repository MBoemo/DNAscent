import sys
import matplotlib
matplotlib.use('Agg')
from matplotlib import pyplot as plt
import numpy as np
import matplotlib.mlab as mlab
from scipy import stats
from scipy import special
import math

#outlier detection
from sklearn.cluster import DBSCAN
from sklearn import mixture

f = open(sys.argv[1],'r')

allScores = []

startingPos = -1
gap = 0
attempts = 0
calls = 0
regionBuffer = []

p = 0.2 #0.3
resolution = 2000
first = True
regionThreshold = -2
printThisOne = False
license = False

for line in f:

	splitLine = line.split('\t')

	if line[0] == '>':

		if not first:

			if len(regionBuffer) >= 10 and printThisOne:

				f = open(readID[1:] + '.bedgraph','w')

				f.write( 'track type=bedGraph name="'+readID+'" description="BedGraph format" visibility=full color=200,100,0 altColor=0,100,200 priority=20'+'\n')
				f.write( " ".join(regionBuffer[0][0:4])+'\n' )

				for i, bufEntry in enumerate(regionBuffer[1:-1]):

					if bufEntry[4] == "BrdU" and regionBuffer[i][4] == "Thym" and regionBuffer[i+2][4] == "Thym":

						bufEntry[4] = "Thym"

					elif bufEntry[4] == 'Thym' and regionBuffer[i][4] == 'BrdU' and regionBuffer[i+2][4] == 'BrdU':

						bufEntry[4] = 'BrdU'

					allScores.append(float(bufEntry[3]))

					f.write( " ".join(bufEntry[0:4])+'\n')
				f.write( " ".join(regionBuffer[-1:][0][0:4]) +'\n')
				f.close()

		splitLine_space = line.split(' ')
		splitLine_colon = splitLine_space[1].split(':')
		thisChromosome = splitLine_colon[0]
		readID = splitLine_space[0]
		startingPos = -1
		gap = 0
		attempts = 0
		calls = 0
		regionBuffer = []
		first = False
		printThisOne = False
		license = False

	else:

		pos = int(splitLine[0])
		pos = pos % 9136
		if abs(pos - 2000) < 21 or abs(pos - 4000) < 21 or abs(pos - 6000) < 21 or abs(pos - 8000) < 21:

			license = True
		
		if license:

			if float(splitLine[1]) > 2.5:
		
				calls += 1

			if startingPos == -1:
				startingPos = int(splitLine[0])

			gap = int(splitLine[0]) - startingPos

			attempts += 1
			if gap > resolution:

				zScore = (calls - attempts*p) / math.sqrt( attempts*p*(1 - p) )

				tempCall = "Thym"

				if zScore >= regionThreshold:
					tempCall = "BrdU"
					printThisOne = True

				regionBuffer.append( [thisChromosome, str(startingPos), splitLine[0], str(zScore+abs(regionThreshold)), tempCall] )
				startingPos = -1
				gap = 0
				attempts = 0
				calls = 0

f.close()

plt.figure(1)
plt.hist(allScores, 35, linewidth=0)
plt.legend(framealpha=0.5)
plt.xlim(-10, 10)
plt.xlabel('Z-Score (Absolute Value)')
plt.ylabel('Count')
plt.title('Distribution of Binomial Z-Scores')
plt.savefig('binomial_zscores.pdf')

