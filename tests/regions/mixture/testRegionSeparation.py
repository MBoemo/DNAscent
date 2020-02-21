import matplotlib
matplotlib.use('Agg')
from matplotlib import pyplot as plt
import numpy as np
import sys

#Usage: python testClustering.py DNAscent.regions
#where DNAscent.regions is the normal output from DNAscent regions

regionScores = []
f = open(sys.argv[1],'r')
for line in f:

	if len(line.rstrip()) == 0:
		continue

	if line[0] == '>':
		continue
	else:
		splitLine = line.rstrip().split()
		regionScores.append(float(splitLine[2]))
f.close()

plt.figure()
plt.hist(regionScores,50)
plt.xlabel('Region Score')
plt.ylabel('Count')
plt.savefig('regionScoreHistogram.pdf')
plt.close()
