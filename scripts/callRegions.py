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

f = open(sys.argv[1],'r')
output = []

target = 'dbscan_ligation_25percentBrdU'

first = True
turnOn = True
for line in f:

	splitLine = line.split('\t')

	if line[0] == '>':

		splitLine_space = line.split(' ')
		splitLine_colon = splitLine_space[1].split(':')
		thisChromosome = splitLine_colon[0]

		#if thisChromosome == target:
		#	turnOn = True
		#else:
		#	turnOn = False


		if not first and len(positions) > 0:
		
			ar = np.array(positions)
			db = DBSCAN( eps = 100, min_samples = 4 ).fit(ar.reshape(-1,1))
			singularities_filtered = []
			for j, label in enumerate(db.labels_):
				if label == -1:
					continue
				else:
					output.append(positions[j])
#					singularities_filtered.append(positions[j])

#			if len(singularities_filtered)>0:
#				ar = np.array(singularities_filtered)
#				db = DBSCAN( eps = 1000, min_samples = 25 ).fit(ar.reshape(-1,1))
#				singularities_filtered = []
#				for j, label in enumerate(db.labels_):
#					if label == -1:
#						continue
#					else: 
#						output.append(positions[j])

			


		positions = []
		first = False

	elif turnOn:

		if float(splitLine[1]) > 2.5:
			positions.append(int(splitLine[0]))

plt.hist(output,400,alpha=0.5,linewidth=0)
plt.title('Ligation Substrate - 25% BrdU')
plt.xlabel('Position on Chromosome (bp)')
plt.ylabel('BrdU Call Count')
plt.savefig(target+'.pdf')
f.close()
