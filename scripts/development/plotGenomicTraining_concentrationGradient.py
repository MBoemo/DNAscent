import sys
import matplotlib
matplotlib.use('Agg')
from matplotlib import pyplot as plt
import numpy as np
import matplotlib.mlab as mlab
from scipy import stats

totalReads = 10000

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

#eventalign output 0%
sixmer2eventsBrdU_0 = {}
f = open(sys.argv[1],'r')
currentRead = ''
readCounter = 0
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

	if readCounter == totalReads:
		break

	eventTime = float(splitLine[8])
	if eventTime < 0.002:
		continue

	sixmer = splitLine[9]
	if sixmer not in sixmer2eventsBrdU_0:
		sixmer2eventsBrdU_0[sixmer] = []
	else:
		sixmer2eventsBrdU_0[sixmer].append( float(splitLine[6]) )
f.close()

#eventalign output 60%
sixmer2eventsBrdU_60 = {}
f = open(sys.argv[2],'r')
currentRead = ''
readCounter = 0
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

	if readCounter == totalReads:
		break

	eventTime = float(splitLine[8])
	if eventTime < 0.002:
		continue

	sixmer = splitLine[9]
	if sixmer not in sixmer2eventsBrdU_60:
		sixmer2eventsBrdU_60[sixmer] = []
	else:
		sixmer2eventsBrdU_60[sixmer].append( float(splitLine[6]) )
f.close()

#eventalign output 80%
sixmer2eventsBrdU_80 = {}
f = open(sys.argv[3],'r')
currentRead = ''
readCounter = 0
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

	if readCounter == totalReads:
		break

	eventTime = float(splitLine[8])
	if eventTime < 0.002:
		continue

	sixmer = splitLine[9]
	if sixmer not in sixmer2eventsBrdU_80:
		sixmer2eventsBrdU_80[sixmer] = []
	else:
		sixmer2eventsBrdU_80[sixmer].append( float(splitLine[6]) )
f.close()

#eventalign output 100%
sixmer2eventsBrdU_100 = {}
f = open(sys.argv[4],'r')
currentRead = ''
readCounter = 0
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

	if readCounter == totalReads:
		break

	eventTime = float(splitLine[8])
	if eventTime < 0.002:
		continue

	sixmer = splitLine[9]
	if sixmer not in sixmer2eventsBrdU_100:
		sixmer2eventsBrdU_100[sixmer] = []
	else:
		sixmer2eventsBrdU_100[sixmer].append( float(splitLine[6]) )
f.close()

#pore model
model = {}
f = open(sys.argv[5],'r')
for line in f:

	if line[0] == '#' or line[0] == 'k':
		continue
	
	splitLine = line.rstrip().split('\t')
	model[splitLine[0]] = [	float(splitLine[1]), float(splitLine[2]) ]
f.close()



for i, key in enumerate(sixmer2eventsBrdU_60):

	if key == 'NNNNNN':
		continue
	
	#if key in mixture:
	if len( sixmer2eventsBrdU_0[key] ) > 100 and len( sixmer2eventsBrdU_60[key] ) > 100 and len( sixmer2eventsBrdU_80[key] ) > 100 and len( sixmer2eventsBrdU_100[key] ) > 100:
		x = np.linspace( np.mean(sixmer2eventsBrdU_60[key])-15, np.mean(sixmer2eventsBrdU_60[key])+15, 1000 )
		plt.figure(i)

		#density
		densityBrdU = stats.kde.gaussian_kde( sixmer2eventsBrdU_0[key] )
		plt.plot(x, densityBrdU(x), label='0% BrdU Density')
		densityBrdU = stats.kde.gaussian_kde( sixmer2eventsBrdU_60[key] )
		plt.plot(x, densityBrdU(x), label='30% BrdU Density')
		densityBrdU = stats.kde.gaussian_kde( sixmer2eventsBrdU_80[key] )
		plt.plot(x, densityBrdU(x), label='50% BrdU Density')
		densityBrdU = stats.kde.gaussian_kde( sixmer2eventsBrdU_100[key] )
		plt.plot(x, densityBrdU(x), label='80% BrdU Density')
		
		#pore model
		yModel = mlab.normpdf( x, model[key][0], model[key][1] )
		plt.plot(x, yModel, label='6mer Pore Model')

		#plotting stuff
		plt.xlabel('pA')
		plt.ylabel('Count')
		plt.title( key )
		plt.legend(loc='upper right')
		plt.savefig( key + '.png' )
		plt.close()
