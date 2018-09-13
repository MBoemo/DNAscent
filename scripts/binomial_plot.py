import matplotlib
matplotlib.use('Agg')
from matplotlib import pyplot as plt
import sys
import numpy as np

f = open(sys.argv[1],'r')
target = 'chrIV'
chrLength = 1531933
regionCounts = np.array([0]*chrLength)

for line in f:

	if line[0] == '#':
		continue

	splitLine = line.rstrip().split(' ')
	
	if splitLine[0] == target:
		if float(splitLine[3]) > -1 and float(splitLine[3]) < 2:
			regionCounts[int(splitLine[1]):int(splitLine[2])] += 1

f.close()

xArr = np.array(range(0,chrLength))

plt.plot(xArr[0::5000], regionCounts[0::5000] )
plt.xlabel('Location on Chromosome (bp)')
plt.ylabel('Count')
plt.title('BrdU Calls in Chromosomal Regions')
plt.savefig('regions_plot.pdf')
