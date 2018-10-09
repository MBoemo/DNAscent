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
chrOrigins2IP = {}
early_chromosome2origins = []
late_chromosome2origins = []
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

		early_chromosome2origins.append( (chromosome, oriLoc) )

	elif splitLine[4] == "0": #checked (late)

		late_chromosome2origins.append( (chromosome, oriLoc) )

	else:
		print "indexing issue"

	#keep track of efficiency at each origin
	efficiency = float(splitLine[5])
	chrOrigins2efficiency[(chromosome,oriLoc)] = efficiency	
	chrOrigins2IP[(chromosome,oriLoc)] = float(splitLine[6])			


f.close()
print "bed imported"


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
print "regions file imported"

#create the heatmap matrix
yEarly = [[] for _ in range(2*heatmapWidth+1)]
yLate = [[] for _ in range(2*heatmapWidth+1)]

matrix = np.ones((rowCount, 2*heatmapWidth+1)) * np.nan
checkedHM = np.ones((rowCount, 100)) * np.nan
efficiencyHM = np.ones((rowCount, 100)) * np.nan
BrdUIPHM = np.ones((rowCount, 100)) * np.nan

row = 0
for t in sorted_chrOrigins2efficiency:

	chromOriPair = t[0]
	efficiency = t[1]

	if chromOriPair in chrOrigins2candidiates:
		for read in chrOrigins2candidiates[chromOriPair]:
			originPos = chromOriPair[1]

			for startEndScore in read:

				if (originPos - startEndScore[0]) < heatmapWidth or (startEndScore[1] - originPos) < heatmapWidth:

					#add to heatmap
					lower = startEndScore[0] - originPos + heatmapWidth
					higher = startEndScore[1] - originPos + heatmapWidth
					if lower > 0 and higher < 2*heatmapWidth + 1:

						matrix[row, lower-100:higher+100] = startEndScore[2] #main heatmap
						efficiencyHM[row, :] = chrOrigins2efficiency[chromOriPair] #efficiency heatmap
						BrdUIPHM[row, :] = chrOrigins2IP[chromOriPair] #ip heatmap

						#add to line plots
						if chromOriPair in late_chromosome2origins:

							checkedHM[row, :] = 1 #checked heatmap

							for i in range(lower,higher):

								yLate[i].append(startEndScore[2])

						elif chromOriPair in early_chromosome2origins:

							checkedHM[row, :] = 0 #unchecked heatmap

							for i in range(lower,higher):

								yEarly[i].append(startEndScore[2])

			row += 1
print "heatmap matrices made"

#statistics for the line plots

yEarly_avg = []
yEarly_stdv = []
xEarly = []
for i in range(0, len(yEarly), 100):
	if len(yEarly[i]) == 0:
		continue
	yEarly_avg.append( np.mean(yEarly[i]) )
	yEarly_stdv.append( np.std(yEarly[i]) )
	xEarly.append(i)
yLate_avg = []
yLate_stdv = []
xLate = []
for i in range(0, len(yLate), 100):
	if len(yLate[i]) == 0:
		continue
	yLate_avg.append( np.mean(yLate[i]) )
	yLate_stdv.append( np.std(yLate[i]) )
	xLate.append(i)

yLate_avg = np.array( yLate_avg )
yLate_stdv = np.array( yLate_stdv )
yEarly_avg = np.array( yEarly_avg )
yEarly_stdv = np.array( yEarly_stdv )
print "line plot statistics computed"

#plot the main heatmap
plt.figure(1)
cmap = plt.cm.jet#plt.cm.copper_r
cmap.set_bad(color='w')
plt.imshow(matrix[0:rowCount,:],cmap=cmap,aspect = 'auto',interpolation='none')
plt.colorbar()
plt.clim(-5,1)
plt.xticks([0, 5000, 10000, 15000, 20000, 25000, 30000],('-15', '-10', '-5', '0', '5', '10', '15'))
plt.xlabel('Distance from Origin (kb)')
plt.ylabel('Individual Reads')
plt.savefig('OriHeatmap_efficiencyRanked.pdf')
plt.clf()
plt.cla()
plt.close()
print "main heatmap plotted"

#efficiency heatmap
plt.figure(2)
cmap = plt.cm.copper
cmap.set_bad(color='w')
plt.imshow(efficiencyHM[0:rowCount,:],cmap=cmap,aspect = 'auto',interpolation='none')
plt.colorbar()
plt.savefig('OriHeatmap_efficiencyHM.pdf')
plt.clf()
plt.cla()
plt.close()
print "efficiency heatmap plotted"

#checked heatmap
plt.figure(3)
cmap = plt.cm.copper
cmap.set_bad(color='w')
plt.imshow(checkedHM[0:rowCount,:],cmap=cmap,aspect = 'auto',interpolation='none')
plt.colorbar()
plt.savefig('OriHeatmap_checkedHM.pdf')
plt.clf()
plt.cla()
plt.close()
print "checked heatmap plotted"

#IP heatmap
plt.figure(4)
cmap = plt.cm.copper
cmap.set_bad(color='w')
plt.imshow(BrdUIPHM[0:rowCount,:],cmap=cmap,aspect = 'auto',interpolation='none')
plt.colorbar()
plt.savefig('OriHeatmap_ipHM.pdf')
plt.clf()
plt.cla()
plt.close()
print "IP heatmap plotted"

#make the line plot
plt.figure(5)
plt.plot( xEarly, yEarly_avg,'k',label='Unchecked')
print len(xEarly), len(yEarly_avg), len(yEarly_stdv)
plt.fill_between(xEarly, yEarly_avg - yEarly_stdv, yEarly_avg + yEarly_stdv, color='k',alpha=0.3,linewidth=0 )
plt.plot( xLate, yLate_avg,'r',label='Checked')
plt.fill_between(xLate, yLate_avg - yLate_stdv, yLate_avg + yLate_stdv, color='r',alpha=0.3,linewidth=0 )
plt.xticks([0, 5000, 10000, 15000, 20000, 25000, 30000],('-15', '-10', '-5', '0', '5', '10', '15'))
plt.legend()
plt.xlabel('Distance from Origin (kb)')
plt.ylabel('Z-Score')
plt.savefig('OriHeatmap_linePlot.pdf')
