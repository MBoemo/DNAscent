import matplotlib
matplotlib.use('Agg')
from matplotlib import pyplot as plt
import sys
import numpy as np

#path to oridb bed file
f_oridb = '/home/mb915/rds/rds-mb915-notbackedup/oridb.bed'

def namefix(num):

	if num == 'chr1':
		return 'chrI'
	elif num == 'chr2':
		return 'chrII'
	elif num == 'chr3':
		return 'chrIII'
	elif num == 'chr4':
		return 'chrIV'
	elif num == 'chr5':
		return 'chrV'
	elif num == 'chr6':
		return 'chrVI'
	elif num == 'chr7':
		return 'chrVII'
	elif num == 'chr8':
		return 'chrVIII'
	elif num == 'chr9':
		return 'chrIX'
	elif num == 'chr10':
		return 'chrX'
	elif num == 'chr11':
		return 'chrXI'
	elif num == 'chr12':
		return 'chrXII'
	elif num == 'chr13':
		return 'chrXIII'
	elif num == 'chr14':
		return 'chrXIV'
	elif num == 'chr15':
		return 'chrXV'
	elif num == 'chr16':
		return 'chrXVI'
	else:
		print(num)
		print('problem')

#parse oridb
f = open(f_oridb,'r')
oridb_chr2bounds = {}
for line in f:
	splitLine = line.rstrip().split()
	chromosome = namefix(splitLine[0])
	lb = int(splitLine[1])
	ub = int(splitLine[2])
	if chromosome in oridb_chr2bounds:
		oridb_chr2bounds[chromosome].append((lb,ub))
	else:
		oridb_chr2bounds[chromosome] = [(lb,ub)]
f.close()

#parse the origin calls
distances = []
f = open(sys.argv[1],'r')
for line in f:
	splitLine = line.rstrip().split()
	chromosome = splitLine[0]

	if chromosome == 'chrM':
		continue

	lb = int(splitLine[2]) 
	ub = int(splitLine[1])	

	minDist = 1000000000
	for ori in oridb_chr2bounds[chromosome]:
		if ori[0] < lb < ori[1] or ori[0] < ub < ori[1] or lb < ori[0] < ub or lb < ori[1] < ub:
			minDist = 0
			break
		elif (ori[1] < lb and lb - ori[1] < minDist):
			minDist = lb - ori[1]
		elif (ub < ori[0] and ori[0] - ub < minDist):
			minDist = ori[0] - ub
	distances.append(minDist)
f.close()	

plt.figure()
plt.hist(distances,50)
plt.xlabel('Distance to Nearest Origin')
plt.ylabel('Count')
plt.savefig('distToNearestOri.pdf')
plt.close()

print(np.mean(distances))
