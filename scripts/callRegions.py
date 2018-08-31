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
f_out = open(sys.argv[1] + '.dbscanFiltered','w')
output = []
lineHolder = []
idx = 0
indices = []

first = True
for line in f:

	splitLine = line.split('\t')

	if line[0] == '>':

		splitLine_space = line.split(' ')
		splitLine_colon = splitLine_space[1].split(':')
		thisChromosome = splitLine_colon[0]

		if not first and len(positions) > 0:
		
			ar = np.array(positions)
			db = DBSCAN( eps = 100, min_samples = 2 ).fit(ar.reshape(-1,1))
			for j, label in enumerate(db.labels_):
				if label == -1:
					lineHolder[indices[j]] = (lineHolder[indices[j]][0], 0, 'f')
				#else:
				#	f_out.write(lineHolder[j])

			for t in lineHolder:
				f_out.write(str(t[0]) + '\t' + str(t[1]) + '\t' + str(t[2]) + '\n')

		positions = []
		lineHolder = []
		indices = []
		idx = 0
		first = False
		f_out.write(line)

	else:

		if float(splitLine[1]) > 2.5:

			positions.append(int(splitLine[0]))
			indices.append(idx)

		if float(splitLine[1]) > 2.5:
			call = 1
		else:
			call = 0
		
		lineHolder.append( (splitLine[0], call, '') )
		idx += 1

f_out.close()
f.close()
