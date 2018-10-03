import sys
import matplotlib
matplotlib.use('Agg')
from matplotlib import pyplot as plt

f = open(sys.argv[1],'r')
positionMap = {0:[],1:[],2:[],3:[],4:[],5:[]}
for line in f:
	splitLine = line.split('\t')

	#if T is in the 6mer	
	if splitLine[0].count('T') == 1:
		positionMap[splitLine[0].find('T')].append(float(splitLine[9]))

plt.figure()
for i in range(0,6):
	plt.subplot(2,3,i+1)
	plt.tick_params(axis='both', which='major', labelsize=4)
	plt.tick_params(axis='both', which='minor', labelsize=4)
	plt.hist(positionMap[i],20,linewidth=0)
	plt.xlim(0,3.5)
	if i in [0,3]:
		plt.ylabel('Count')
	if i in [3,4,5]:
		plt.xlabel('Expected Log Likelihood')


	plt.title('Position '+str(i+1))

plt.tight_layout()
plt.savefig('positionEffects.pdf')
