#----------------------------------------------------------
# Copyright 2017 University of Oxford
# Written by Michael A. Boemo (michael.boemo@path.ox.ac.uk)
# This software is licensed under GPL-3.0.  You should have
# received a copy of the license with this software.  If
# not, please Email the author.
#----------------------------------------------------------

#development script
#measures thymidine training against the ONT pore model

import sys
import matplotlib
matplotlib.use('Agg')
from matplotlib import pyplot as plt
import numpy as np


#MAIN--------------------------------------------------------------------------------------------------------------------------------------
f = open( sys.argv[1],'r' )
g = f.readlines()
f.close()

diffs = []
for line in g:
	if line[0] != '#' and line[0:4] != 'kmer': #ignore the header
		splitLine = line.split('\t')
		if float(splitLine[2]) < 3.5:
			diffs.append( float(splitLine[1]) - float(splitLine[3]) )
g = None

plt.hist(diffs,len(diffs)/10)
plt.xlabel('Difference Between Means (pA)')
plt.ylabel('Count')
plt.title('Comparison Between Trained Thymidine Mean and ONT Mean, N=' + str(len(diffs)))
plt.savefig('modelDifferencePlot.pdf')
print 'Mean: ',np.mean(diffs)
print 'Stdv: ',np.std(diffs)
