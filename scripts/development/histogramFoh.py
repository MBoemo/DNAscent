import sys
import matplotlib
matplotlib.use('Agg')
from matplotlib import pyplot as plt

inputStr = sys.argv[1]

target7mer = inputStr[inputStr.find('_')+2:inputStr.find('.')-1]

f = open(inputStr,'r')
g = f.readlines()
f.close()

sevenMers = []
readCounts = []

first = True 

for line in g:
	if line[0] == '>':

		if not first:
			sevenMers.append(current7mer)
			readCounts.append(count)

		current7mer = line[1:].rstrip()
		count = 0
		first = False		

	elif line[0] in ['A','T','G','C']:
		count += 1

sevenMers.append(current7mer)
readCounts.append(count)

sorted_readCounts = sorted(readCounts)
sorted_sevenMers = [x for _,x in sorted(zip(readCounts,sevenMers))]

sorted_readCounts = sorted_readCounts[::-1]
sorted_sevenMers = sorted_sevenMers[::-1]

upTo = min([len(sorted_readCounts),20]) #plot the top 20

plt.bar(range(0,upTo),sorted_readCounts[0:upTo],align='edge')
plt.xticks(range(0,upTo), sorted_sevenMers[0:upTo], rotation='vertical',fontsize=6)
plt.title('Target Sequence: '+target7mer)

plt.savefig(inputStr+'_barGraph.pdf')


