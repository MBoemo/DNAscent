import sys
import matplotlib
matplotlib.use('Agg')
from matplotlib import pyplot as plt
import numpy as np
import matplotlib.mlab as mlab
from scipy import stats
import math

#outlier detection
from sklearn.cluster import DBSCAN
from sklearn import mixture

target = 'chrVI'
target2 = 'chr6'

#import the BrdU IP data
xIP = []
yIP = []
f = open(sys.argv[2],'r')
for line in f:
	splitLine = line.rstrip().split('\t')
	if splitLine[0] == target2:
		xIP.append( (float(splitLine[2]) + float(splitLine[1]))/2 )
		yIP.append( float(splitLine[3]) )
f.close()

#import the timing data
xT = []
yT = []
f = open(sys.argv[3],'r')
for line in f:
	splitLine = line.rstrip().split('\t')
	if splitLine[0] == target2:
		xT.append( (float(splitLine[2]) + float(splitLine[1]))/2 )
		yT.append( float(splitLine[3]) )
f.close()




f = open(sys.argv[1],'r')
output = []

first = True
for line in f:

	splitLine = line.split('\t')

	if line[0] == '>':

		splitLine_space = line.split(' ')
		splitLine_colon = splitLine_space[1].split(':')
		thisChromosome = splitLine_colon[0]

		if thisChromosome == target:
			turnOn = True
		else:
			turnOn = False
	
		if not first and len(positions) > 0:
		
			ar = np.array(positions)
			db = DBSCAN( eps = 100, min_samples = 4 ).fit(ar.reshape(-1,1))
			singularities_filtered = []
			for j, label in enumerate(db.labels_):
				if label == -1:
					continue
				else: 
					singularities_filtered.append(positions[j])

			if len(singularities_filtered)>0:
				ar = np.array(singularities_filtered)
				db = DBSCAN( eps = 1000, min_samples = 25 ).fit(ar.reshape(-1,1))
				singularities_filtered = []
				for j, label in enumerate(db.labels_):
					if label == -1:
						continue
					else: 
						output.append(positions[j])

			


		positions = []
		first = False

	elif turnOn:

		if float(splitLine[1]) > 2.5:
			positions.append(int(splitLine[0]))

#plt.hist(output,400,alpha=0.5,linewidth=0)
plt.figure()
height, edges = np.histogram(output, bins=400)
histMean = np.mean(height)
ipMean = np.mean(yIP)
tMean = np.mean(yT)
yIP = ( np.array(yIP) ) * (histMean/ipMean)
yT = ( np.array(yT) ) * (float(max(height))/float(max(yT)))
center = (edges[:-1] + edges[1:]) / 2
#plt.bar(center, height, align='center', alpha=0.5, linewidth=0)
print center
print height
plt.bar(center, height, align='center',color='k', alpha=0.5,label='Nanopore BrdU Calls')
plt.plot(xIP,yIP,label='BrdU IP')
plt.plot(xT,yT,label='Replication Timing')
plt.legend()
plt.title(target)
plt.xlabel('Position on Chromosome (bp)')
plt.ylabel('BrdU Call Count')
plt.savefig(target+'.pdf')
f.close()
