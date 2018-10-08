import sys
import numpy as np
import matplotlib
matplotlib.use('Agg')
from matplotlib import pyplot as plt

target = 'chrII'
length = 813184
BrdUCalls = [0]*length
coverage = [0]*length
BrdUCalls = np.array(BrdUCalls)
coverage = np.array(coverage)

#import IP data
xIP = []
yIP = []
f = open(sys.argv[1], 'r')
turnOn = False

for line in f:

	splitLine = line.rstrip().split('\t')
	thisChromosome = splitLine[0]

	if thisChromosome == target:

		xIP.append((float(splitLine[1]) + float(splitLine[2])) / 2.0 )
		yIP.append(float(splitLine[4]))
f.close()
yIP = np.convolve(yIP, np.ones((10,))/10, mode='same')


#import detect data
f = open(sys.argv[2],'r')
for line in f:

	if line[0] == '>':

		splitLine = line.rstrip().split(' ')
		splitLine2 = splitLine[1].split(':')
		chromosome = splitLine2[0]

		continue

	else:

		splitLine = line.rstrip().split('\t')

		if chromosome == target:

			coverage[int(splitLine[0])] += 1
			if float(splitLine[1]) > 2.5:

				BrdUCalls[int(splitLine[0])] += 1



xBrdU = []
yBrdU = []
for i in range( 0, length, 100 ):

	if float(sum( coverage[i:i+100])) == 0.0:
		continue
	else:
		yBrdU.append(float(sum( BrdUCalls[i:i+100] )) / float(sum( coverage[i:i+100])))
		xBrdU.append(i+50)

yBrdUSmooth = np.convolve(yBrdU, np.ones((10,))/10, mode='same')


f = open('BrdU_data.plot','w')
for i, v in enumerate(yBrdUSmooth):
	f.write(str(xBrdU[i]) + ' ' + str(v) + '\n')
f.close()
		

#normalise for axes scales
yBrdUSmooth = np.array(yBrdUSmooth)
normFactor = np.mean(yIP) / np.mean(yBrdUSmooth)
yBrdUSmooth = yBrdUSmooth*normFactor


plt.figure(1)
plt.plot(xIP, yIP, label='IP', alpha=0.5)
plt.plot(xBrdU, yBrdUSmooth, label='Nanopore', alpha=0.5)
plt.xlim(0,length)
#plt.ylim(0,6)
plt.legend(framealpha=0.5)
plt.xlabel('Position on Chromosome (bp)')
plt.ylabel('A.U.')
plt.title(target + ': BrdU IP vs. Nanopore Region Calls')
plt.savefig('IP_vs_nanopore.pdf')

plt.figure(2)
plt.plot(range(0,length), coverage)
plt.xlabel('Position on Chromosome (bp)')
plt.ylabel('Count')
plt.title(target + ': Coverage')
plt.savefig('coverage.pdf')
