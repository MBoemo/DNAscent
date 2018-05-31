import sys
import matplotlib
matplotlib.use('Agg')
from matplotlib import pyplot as plt
import numpy as np
import matplotlib.mlab as mlab
from scipy import stats



#first argument is the reference, second argument is the thymidine working data, and third argument is the 6mer pore model
#--------------------------------------------------------------------------------------------------------------------------------------
def import_reference(filename):

	f = open(filename,'r')
	g = f.readlines()
	f.close()

	reference = ''
	for line in g:
		if line[0] != '>':
			reference += line.rstrip()
	g = None

	reference = reference.upper()

	if not all(c in ['A','T','G','C'] for c in reference):
		warnings.warn('Warning: Illegal character in reference.  Legal characters are A, T, G, and C.', Warning)

	return reference

#--------------------------------------------------------------------------------------------------------------------------------------
reference = import_reference( sys.argv[1] )

f = open(sys.argv[2],'r')
pos2eventsT = {}
for line in f:

	splitLine = line.rstrip().split(' ')
	position = splitLine[0]
	if position not in pos2eventsT:
		pos2eventsT[position] = []
	else:
		pos2eventsT[position] += map(float,splitLine[1:])
f.close()

#--------------------------------------------------------------------------------------------------------------------------------------
model = {}
f = open(sys.argv[3],'r')
for line in f:

	if line[0] == '#' or line[0] == 'k':
		continue
	
	splitLine = line.rstrip().split('\t')
	model[splitLine[0]] = [	float(splitLine[1]), float(splitLine[2]) ]
f.close()

#--------------------------------------------------------------------------------------------------------------------------------------
for i, key in enumerate(pos2eventsT):
	
	sixMer = reference[int(key):int(key)+6]
	x = np.linspace( np.mean(pos2eventsT[key])-15, np.mean(pos2eventsT[key])+15, 1000 )
	plt.figure(i)
	densityT = stats.kde.gaussian_kde( pos2eventsT[key] )
	plt.plot( x, densityT(x), label='Thymidine Density')
	yModel = mlab.normpdf( x, model[sixMer][0], model[sixMer][1] )
	plt.plot( x, yModel, label='6mer Pore Model')
	plt.xlabel('pA')
	plt.ylabel('Count')
	plt.title( reference[int(key):int(key)+6] )
	plt.legend(loc='upper right')
	plt.savefig( sixMer + '.png' )
	plt.close()
