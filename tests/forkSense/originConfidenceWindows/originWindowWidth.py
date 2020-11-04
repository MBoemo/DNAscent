import matplotlib
matplotlib.use('Agg')
from matplotlib import pyplot as plt
import sys

buf = []

f = open(sys.argv[1],'r')
for line in f:
	splitLine = line.rstrip().split()
	chromosome = splitLine[0]
	if chromosome == 'chrM':
		continue
	dist = int(splitLine[2]) - int(splitLine[1])
	buf.append(dist)

f.close()

plt.figure()
plt.hist(buf,50)
plt.xlim(0,10000)
plt.xlabel('Size of Origin Call Confidence Window')
plt.ylabel('Count')
plt.savefig('originWidths.pdf')
plt.close()
