import sys
import matplotlib
matplotlib.use('Agg')
from matplotlib import pyplot as plt
import numpy as np

repeatLen = 9136
coverage = np.array([0]*repeatLen)
scores = np.array([0.0]*repeatLen)

f = open(sys.argv[1],'r')
for line in f:

	if line[0] == '>':
		continue

	else:
		splitLine = line.rstrip().split('\t')

		lower = int(splitLine[0]) % repeatLen
		upper = int(splitLine[1]) % repeatLen
		score = float(splitLine[2])

		if upper < lower:
			if score > 0:
				coverage[lower:repeatLen] += 1
				scores[lower:repeatLen] += score
				coverage[0:upper] += 1
				scores[0:upper] += score
		else:
			if score > 0:
				coverage[lower:upper] += 1
				scores[lower:upper] += score

f.close()

outMat = np.divide(scores,coverage,where=coverage>0)

plt.figure()
plt.plot(range(0,repeatLen),np.concatenate((outMat[2932:], outMat[0:2932]),axis=None))
plt.scatter([4568], [outMat[2932]],c='r')
#plt.ylim(-1, 0)
plt.xlim(0,9136)
plt.xlabel('Position on rDNA Repeat (bp)')
plt.ylabel('Score')
plt.savefig('rDNA.pdf')
plt.close()

plt.clf()
plt.figure()
plt.plot(range(0,repeatLen),coverage)	
plt.xlim(0,9136)
plt.xlabel('Position on rDNA Repeat (bp)')
plt.ylabel('Coverage')
plt.savefig('coverage.pdf')
plt.close()
