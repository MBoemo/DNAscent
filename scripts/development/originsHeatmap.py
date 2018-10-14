import sys
import os
import numpy as np
import matplotlib
matplotlib.use('Agg')
from matplotlib import pyplot as plt
import operator
import seaborn as sns
from mpl_toolkits.mplot3d import Axes3D

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
sorted_chrOrigins2IP = sorted(chrOrigins2IP.items(), key=operator.itemgetter(1))
sorted_chrOrigins2IP = sorted_chrOrigins2IP[::-1] #sort from high to low efficiency

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
		splitLine4=splitLine3[1].split('#')
		readEnd = int(splitLine4[0])
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
			if oriPos > lower and oriPos < upper and call == "BrdU" and abs(readStart - oriPos) > 4000 and abs(readEnd - oriPos) > 4000:

				centeredOrigins.append(oriPos)
print "regions file imported"

#create the heatmap matrix
yEarly = [[] for _ in range(2*heatmapWidth+1)]
yLate = [[] for _ in range(2*heatmapWidth+1)]

matrix = np.ones((rowCount, 2*heatmapWidth+1)) * np.nan
checkedHM = np.ones((rowCount, 100)) * np.nan
efficiencyHM = np.ones((rowCount, 100)) * np.nan
BrdUIPHM = np.ones((rowCount, 100)) * np.nan

symmetry_left = []
symmetry_right = []

row = 0
for t in sorted_chrOrigins2IP:

	chromOriPair = t[0]
	efficiency = t[1]

	if chromOriPair in chrOrigins2candidiates:
		for read in chrOrigins2candidiates[chromOriPair]:

			originPos = chromOriPair[1]
			centeredIndex = 0

			for sec, startEndScore in enumerate(read):

				if originPos > startEndScore[0] and originPos < startEndScore[1]:
					centeredIndex = sec
					

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

			#for symmetry heatmap, only do if there are at least 8kb of read on either side of the origin
			if centeredIndex >=8 and len(read) - centeredIndex >= 8:

				l = 0
				r = 0
				for right_startEndScore in read[centeredIndex+1:]:
					score = right_startEndScore[2]
					if score < 0 :
						r = abs(right_startEndScore[1] - read[centeredIndex][1])
						break
				for left_startEndScore in read[0:centeredIndex][::-1]:
					score = left_startEndScore[2]
					if score < 0 :
						l = abs(left_startEndScore[0] - read[centeredIndex][0])
						break
				if l > 0 and r > 0:
					print row
					symmetry_left.append(max(l,r))
					symmetry_right.append(min(l,r))

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
	yEarly_stdv.append( np.std(yEarly[i]) / np.sqrt(len(yEarly[i]) ) )
	xEarly.append(i)
yLate_avg = []
yLate_stdv = []
xLate = []
for i in range(0, len(yLate), 100):
	if len(yLate[i]) == 0:
		continue
	yLate_avg.append( np.mean(yLate[i]) )
	yLate_stdv.append( np.std(yLate[i]) / np.sqrt(len(yLate[i]) ) )
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
plt.clim(-4,4)
plt.xticks([0, 5000, 10000, 15000, 20000, 25000, 30000],('-15', '-10', '-5', '0', '5', '10', '15'))
plt.xlabel('Distance from Origin (kb)')
plt.ylabel('Individual Reads')
plt.savefig('OriHeatmap_readsRankedbyBrdUIP.pdf')
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
plt.fill_between(xEarly, yEarly_avg - yEarly_stdv, yEarly_avg + yEarly_stdv, color='k',alpha=0.3,linewidth=0 )
plt.plot( xLate, yLate_avg,'r',label='Checked')
plt.fill_between(xLate, yLate_avg - yLate_stdv, yLate_avg + yLate_stdv, color='r',alpha=0.3,linewidth=0 )
plt.xticks([0, 5000, 10000, 15000, 20000, 25000, 30000],('-15', '-10', '-5', '0', '5', '10', '15'))
plt.legend()
plt.xlabel('Distance from Origin (kb)')
plt.ylabel('Z-Score')
plt.savefig('OriHeatmap_linePlot.pdf')
plt.clf()
plt.cla()
plt.close()

#make the symmetry scatter plot
sns.kdeplot(symmetry_left, symmetry_right, cmap="Reds", shade=True, bw=1000,gridsize=100)

plt.xlabel('Leftward Distance (bp)')
plt.ylabel('Rightward Distance (bp)')
plt.savefig("OriHeatmap_density.pdf")
plt.clf()
plt.cla()
plt.close()

#scatter plot for symmetry
fig = plt.figure(6)
ax = fig.add_subplot(111, projection='3d')
symmetry_left = np.array(symmetry_left) / 1000
symmetry_right = np.array(symmetry_right) / 1000
hist, xedges, yedges = np.histogram2d(symmetry_left, symmetry_right, bins=8, range=[[0, 16], [0, 16]])
xpos, ypos = np.meshgrid(xedges[:-1] + 0.25, yedges[:-1] + 0.25)
xpos = xpos.flatten('F')
ypos = ypos.flatten('F')
zpos = np.zeros_like(xpos)

# Construct arrays with the dimensions for the 16 bars.
dx = np.ones_like(zpos)
dy = dx.copy()
dz = hist.flatten()

#get rid of zeros
xpos_f = []
ypos_f = []
zpos_f = []
dx_f = []
dy_f = []
dz_f = []
for i, height in enumerate(dz):
	
	if height > 0:
		xpos_f.append(xpos[i])
		ypos_f.append(ypos[i])		 
		zpos_f.append(zpos[i])
		dx_f.append(dx[i])
		dy_f.append(dy[i])
		dz_f.append(dz[i])

ax.bar3d(xpos_f, ypos_f, zpos_f, dx_f, dy_f, dz_f, color='b', zsort='average')
plt.xlabel('Minimum Distance Travelled (kb)')
plt.ylabel('Maximum Distance Travelled (kb)')
ax.view_init(30, 45)
plt.savefig('OriHeatmap_scatter.pdf')
