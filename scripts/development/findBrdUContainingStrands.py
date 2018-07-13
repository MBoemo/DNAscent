import sys
import matplotlib
matplotlib.use('Agg')
from matplotlib import pyplot as plt
import numpy as np
import matplotlib.mlab as mlab
from scipy import stats

totalReads = 100000

#--------------------------------------------------------------------------------------------------------------------------------------
def displayProgress(current, total):

	barWidth = 70
	progress = float(current)/float(total)

	if progress <= 1.0:
		sys.stdout.write('[')
		pos = int(barWidth*progress)
		for i in range(barWidth):
			if i < pos:
				sys.stdout.write('=')
			elif i == pos:
				sys.stdout.write('>')
			else:
				sys.stdout.write(' ')
		sys.stdout.write('] '+str(int(progress*100))+' %\r')
		sys.stdout.flush()
#--------------------------------------------------------------------------------------------------------------------------------------


#eventalign output
sixmer2eventsBrdU = {}
f = open(sys.argv[1],'r')
currentRead = ''
readCounter = 0
linesBuffer = []
scores = []
for line in f:

	splitLine = line.rstrip().split('\t')

	#ignore the header line
	if splitLine[0] == 'contig':
		continue

	readIndex = splitLine[3]
	if readIndex != currentRead:

		currentRead = readIndex
		readCounter += 1
		displayProgress(readCounter, totalReads)

		if len(linesBuffer) > 10:
			diff = 0
			total = 0
			for l in linesBuffer:

				if l[9].find('T') != -1:
					diff += abs(float(l[10]) - float(l[6]))
					total += 1

			scores.append( float(diff)/float(total) )
		linesBuffer = []
	else:
		linesBuffer.append( splitLine )

	if readCounter == totalReads:
		break
f.close()

plt.figure()
plt.hist(scores,100)
plt.xlabel('Mean Difference Between Model and Event')
plt.ylabel('Count')
plt.title('Clustering BrdU-Containing Reads')
plt.savefig('brdu_read_clustering.png')
plt.close()
