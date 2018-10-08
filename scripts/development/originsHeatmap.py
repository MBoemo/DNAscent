import sys
import os
import numpy as np
import matplotlib
matplotlib.use('Agg')
from matplotlib import pyplot as plt
import operator

heatmapWidth = 15000

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


#import origin bed
f = open(sys.argv[1],'r')
originLocations = {}
chrOrigins2efficiency = {}
chrOrigins2candidiates = {}
early_chromosome2origins = {}
late_chromosome2origins = {}
rowCount = 0

#iterate on the bed file
for line in f:

	if line[0] == '#':
		continue

	splitLine = line.rstrip().split('\t')
	chromosome = splitLine[0]
	chromosome = 'chr' + str( int2RomanNumeral(int(chromosome[3:])) ) #convert to roman numerals for SacCer3

	start = int(splitLine[1])
	end = int(splitLine[2])
	oriLoc = (start + end)/2
	if chromosome not in originLocations:
		originLocations[chromosome] = []
	originLocations[chromosome].append(oriLoc)

	#keep track of early or late
	if splitLine[4] == "1": #unchecked (early)

		if chromosome in early_chromosome2origins:

			early_chromosome2origins[chromosome].append( oriLoc )
		else:
		
			early_chromosome2origins[chromosome] = [ oriLoc ]

	elif splitLine[4] == "0": #checked (late)

		if chromosome in late_chromosome2origins:

			late_chromosome2origins[chromosome].append( oriLoc )
		else:
		
			late_chromosome2origins[chromosome] = [ oriLoc ]

	else:
		print "indexing issue"

	#keep track of efficiency at each origin
	efficiency = float(splitLine[5])
	chrOrigins2efficiency[(chromosome,oriLoc)] = efficiency		


f.close()


#sort by efficiency 
sorted_chrOrigins2efficiency = sorted(chrOrigins2efficiency.items(), key=operator.itemgetter(1))
sorted_chrOrigins2efficiency = sorted_chrOrigins2efficiency[::-1] #sort from high to low efficiency

#iterate on the regions file
f = open(sys.argv[2],'r')
first = True
for line in f:

	if line[0] == '>':

		if not first:
		
			for c in centeredOrigins:

				if (chromosome,c) not in chrOrigins2candidiates:

					chrOrigins2candidiates[(chromosome,c)] = []

				chrOrigins2candidiates[(chromosome,c)].append(startEndScore)
				rowCount += 1
				
			

		splitLine = line.rstrip().split(' ')
		splitLine2 = splitLine[1].split(':')
		chromosome = splitLine2[0]
		splitLine3 = splitLine2[1].split('-')
		readStart = int(splitLine3[0])
		readEnd = int(splitLine3[1])
		startEndScore = []
		centeredOrigins = []
		first = False

	else:

		splitLine = line.rstrip().split('\t')
		lower = int(splitLine[0])
		upper = int(splitLine[1])
		score = float(splitLine[2])
		call = splitLine[3]

		startEndScore.append( (lower, upper, score) )

		for oriPos in originLocations[chromosome]:

			#if this entry is a BrdU call that's centered on the origin and the read is of sufficient length, record that origin
			if oriPos > lower and oriPos < upper and call == "BrdU" and abs(readStart - oriPos) > 8000 and abs(readEnd - oriPos) > 8000:

				centeredOrigins.append(oriPos)

#plot the heatmap
matrix = np.ones((rowCount, 2*heatmapWidth+1)) * np.nan
row = 0
for t in sorted_chrOrigins2efficiency:

	chromOriPair = t[0]
	efficiency = t[1]

	if chromOriPair in chrOrigins2candidiates:
		for read in chrOrigins2candidiates[chromOriPair]:
			originPos = chromOriPair[1]

			for startEndScore in read:

				if (originPos - startEndScore[0]) < heatmapWidth or (startEndScore[1] - originPos) < heatmapWidth:
					lower = startEndScore[0] - originPos + heatmapWidth
					higher = startEndScore[1] - originPos + heatmapWidth
					matrix[row, lower-100:higher+100] = startEndScore[2]
			row += 1
	
plt.figure(1)
cmap = plt.cm.jet#plt.cm.copper_r
cmap.set_bad(color='w')
plt.imshow(matrix[0:rowCount,:],cmap=cmap,aspect = 'auto',interpolation='none')
plt.colorbar()
plt.clim(-5,3)
plt.xticks([0, 5000, 10000, 15000, 20000, 25000, 30000],('-15', '-10', '-5', '0', '5', '10', '15'))
plt.xlabel('Distance from Origin (kb)')
plt.ylabel('Individual Reads')
plt.savefig('OriHeatmap_efficiencyRanked.pdf')
plt.clf()
plt.cla()
plt.close()
