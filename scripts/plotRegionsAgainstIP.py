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
	if line[0] == 'v':
		splitLine = line.rstrip().split(' ')
		splitLine2 = splitLine[1].split('=')
		thisChromosome = splitLine2[1]
		if thisChromosome == target:
			turnOn = True
		else:
			turnOn = False
	elif turnOn:
		splitLine = line.rstrip().split('\t')
		xIP.append(float(splitLine[0]))
		yIP.append(float(splitLine[1]))
f.close()

#import detect data
f = open(sys.argv[2],'r')
for line in f:

	if line[0] == '>':
		continue

	splitLine = line.rstrip().split(' ')
	if splitLine[0] == target:

		coverage[int(splitLine[1]):int(splitLine[2])] += 1
		if splitLine[4] == "BrdU":
			BrdUCalls[int(splitLine[1]):int(splitLine[2])] += 1

#normalise for coverage
BrdUCalls = np.divide(BrdUCalls.astype(float), coverage.astype(float), where=coverage != 0)

#normalise for axes scales
normFactor = np.mean(yIP) / np.mean(BrdUCalls)
BrdUCalls = BrdUCalls*normFactor

plt.figure(1)
plt.plot(xIP, yIP, label='IP', alpha=0.5)
plt.plot(range(0,length), BrdUCalls, label='Nanopore', alpha=0.5)
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
