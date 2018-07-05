import sys
import matplotlib
matplotlib.use('Agg')
from matplotlib import pyplot as plt
import numpy as np
import matplotlib.mlab as mlab
from scipy import stats

f = open(sys.argv[1],'r')
likelihood = []
for line in f:
	splitLine = line.split('\t')
	ll = float(splitLine[13])
	if ll != -1:
		likelihood.append(ll)
f.close()
plt.figure()
x = np.linspace( min(likelihood), max(likelihood), 1000 )
densitylikelihood = stats.kde.gaussian_kde( likelihood )
plt.plot(x, densitylikelihood(x), label='BrdU Density')
plt.xlabel('Log Likelihood')
plt.ylabel('Density')
plt.title(sys.argv[1])
plt.savefig( sys.argv[1] + '.png')


