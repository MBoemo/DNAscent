import sys
import os
import numpy as np
import matplotlib
matplotlib.use('Agg')
from matplotlib import pyplot as plt

bandwidth = 15000

#import origin bed
f = open(sys.argv[1],'r')
early = {} #maps chromosome name to origin locations
late = {}
for line in f:
	splitLine = line.rstrip().split('\t')
	
	chromosome = splitLine[0]
	start = int(splitLine[1])
	end = int(splitLine[2])

	if splitLine[4] == "1": #unchecked (early)

		if chromosome in early:

			early[chromosome].append( (start + end)/2 )
		else:
		
			early[chromosome] = [ (start + end)/2 ]

	elif splitLine[4] == "0": #checked (late)

		if chromosome in late:

			late[chromosome].append( (start + end)/2 )
		else:
		
			late[chromosome] = [ (start + end)/2 ]

	else:
		print "indexing issue"
f.close()


#iterate on the bedgraph files
contents = os.listdir(sys.argv[2])
early_chrOrigin_to_scores = {}
late_chrOrigin_to_scores = {}
earlyRowCount = 0
lateRowCount = 0
lateScoreTot = 0
earlyScoreTot = 0
lateScoreCount = 0
earlyScoreCount = 0

for filename in contents:
	
	splitOnDot = filename.split('.')

	if splitOnDot[-1:][0] != "bedgraph":
		continue

	f = open(sys.argv[2] + '/' + filename,'r')
	g = f.readlines()

	BrdUCallCount = 0
	startEndScore = []
	earlyCandidates = []
	lateCandidates = []
	chromosome = ""

	readStart = int(g[1].rstrip().split(' ')[1])
	readEnd = int(g[-1:][0].rstrip().split(' ')[2])		

	for line in g:

		if line[0:5] == "track":
			continue
		
		splitLine = line.rstrip().split(' ')
		
		chromosome = splitLine[0]
		lower = int(splitLine[1])
		upper = int(splitLine[2])
		score = float(splitLine[3])
		call = splitLine[4]

		startEndScore.append( (lower, upper, score) )
		if call == "BrdU":
			BrdUCallCount += 1

		if chromosome in early:
			earlyOriginLocOnChr = early[chromosome]
			for oriPos in earlyOriginLocOnChr:

				if oriPos > lower and oriPos < upper and call == "BrdU" and abs(readStart - oriPos) > 8000 and abs(readEnd - oriPos) > 8000:
					earlyScoreTot += score
					earlyScoreCount += 1
					earlyCandidates.append(oriPos)

		if chromosome in late:
			lateOriginLocOnChr = late[chromosome]
			for oriPos in lateOriginLocOnChr:
				if oriPos > lower and oriPos < upper and call == "BrdU" and abs(readStart - oriPos) > 8000 and abs(readEnd - oriPos) > 8000:
					lateScoreTot += score
					lateScoreCount += 1
					lateCandidates.append(oriPos)
	
	#read must be at least 10kb in length
	if len(startEndScore) < 5:
		continue


	for i in earlyCandidates:
		if BrdUCallCount > 2:
			earlyRowCount += 1
			if (chromosome,i) in early_chrOrigin_to_scores:
				early_chrOrigin_to_scores[(chromosome,i)].append(startEndScore)
			else:
				early_chrOrigin_to_scores[(chromosome,i)] = [startEndScore]

	for i in lateCandidates:
		if BrdUCallCount > 2:
			lateRowCount += 1
			if (chromosome,i) in late_chrOrigin_to_scores:
				late_chrOrigin_to_scores[(chromosome,i)].append(startEndScore)
			else:
				late_chrOrigin_to_scores[(chromosome,i)] = [startEndScore]		

	f.close()

print "early: ", float(earlyScoreTot)/float(earlyScoreCount)
print "late: ", float(lateScoreTot)/float(lateScoreCount)

#plot early
matrix = np.ones((earlyRowCount, 2*bandwidth+1)) * np.nan
row = 0
for key in early_chrOrigin_to_scores:
	chromosome = key[0]
	originPos = key[1]
	
	for read in early_chrOrigin_to_scores[key]:
		for startEndScore in read:
			if (originPos - startEndScore[0]) < bandwidth or (startEndScore[1] - originPos) < bandwidth:
				lower = startEndScore[0] - originPos + bandwidth
				higher = startEndScore[1] - originPos + bandwidth
				matrix[row, lower-100:higher+100] = startEndScore[2]

	row += 1
	
plt.figure(1)
cmap = plt.cm.copper_r
cmap.set_bad(color='w')
plt.imshow(matrix[0:earlyRowCount,:],cmap=cmap,aspect = 'auto',interpolation='none')
plt.colorbar()
plt.clim(-5,3)
plt.xticks([0, 5000, 10000, 15000, 20000, 25000, 30000],('-15', '-10', '-5', '0', '5', '10', '15'))
plt.xlabel('Distance from Origin (kb)')
plt.ylabel('Individual Reads')
plt.savefig('OriHeatmap_early.pdf')
plt.clf()
plt.cla()
plt.close()

#plot late
matrix = np.ones((lateRowCount, 2*bandwidth + 1)) * np.nan
row = 0
for key in late_chrOrigin_to_scores:
	chromosome = key[0]
	originPos = key[1]
	
	for read in late_chrOrigin_to_scores[key]:
		for startEndScore in read:
			if (originPos - startEndScore[0]) < bandwidth or (startEndScore[1] - originPos) < bandwidth:
				lower = startEndScore[0] - originPos + bandwidth
				higher = startEndScore[1] - originPos + bandwidth
				matrix[row, lower-100:higher+100] = startEndScore[2]

	row += 1
	
plt.figure(2)
cmap = plt.cm.copper_r
cmap.set_bad(color='w')
plt.imshow(matrix[0:lateRowCount,:],cmap=cmap,aspect = 'auto',interpolation='none')
plt.colorbar()
plt.clim(-5,3)
plt.xticks([0, 5000, 10000, 15000, 20000, 25000, 30000],('-15', '-10', '-5', '0', '5', '10', '15'))
plt.xlabel('Distance from Origin (kb)')
plt.ylabel('Individual Reads')
plt.savefig('OriHeatmap_late.pdf')

		
