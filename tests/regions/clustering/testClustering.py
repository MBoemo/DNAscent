import matplotlib
matplotlib.use('Agg')
from matplotlib import pyplot as plt
import numpy as np
import sys

#Usage: python testClustering.py DNAscentregions.stderr
#where DNAscentregions.stderr is from stderr after running DNAscent regions with #define TEST_CLUSTERING 1


centroids = []
clusteredPoints = [[],[]]
index = -1

f = open(sys.argv[1],'r')
for line in f:
	if line[0] == '>':
		centroids.append(float(line.rstrip()[1:]))
		index += 1	
	else:
		clusteredPoints[index].append(float(line.rstrip()))
f.close()

plt.figure()
plt.hist(clusteredPoints[1]+clusteredPoints[0],50,alpha=0.1)
plt.hist(clusteredPoints[0],50,alpha=0.3)
plt.scatter(centroids[0],[1])
plt.hist(clusteredPoints[1],50,alpha=0.3)
plt.scatter(centroids[1],[1])
plt.xlabel('P(Thym is BrdU)')
plt.ylabel('Count')
plt.savefig('regions_clusterTest.pdf')
