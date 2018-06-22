import sys
import matplotlib
matplotlib.use('Agg')
from matplotlib import pyplot as plt
import numpy as np
import matplotlib.mlab as mlab
from scipy import stats

totalReads = 100000

#--------------------------------------------------------------------------------------------------------------------------------------
def reverseComplement(sequence):

	#seq to upper case
	sequence = sequence.upper()

	#build the complement
	revComp = ''
	for char in sequence:
		if char == 'A':
			revComp += 'T'
		elif char == 'T':
			revComp += 'A'
		elif char == 'C':
			revComp += 'G'
		elif char == 'G':
			revComp += 'C'
		#guard against illegal characters
		else:
			warnings.warn('Warning: Illegal character in sequence.  Legal characters are A, T, G, and C.', Warning)

	#take the reverse of the complement and return it
	return revComp[::-1]

#--------------------------------------------------------------------------------------------------------------------------------------
def importReference(filename):

	f = open(filename,'r')
	g = f.readlines()
	f.close()

	reference = {}
	for line in g:

		if line[0] == '>':
			currentChromosome = line[1:].rstrip()
			reference[ currentChromosome ] = ''

		if line[0] != '>':
			reference[ currentChromosome ] += line.rstrip().upper()

	return reference

#--------------------------------------------------------------------------------------------------------------------------------------
def displayProgress(current, total):

	barWidth = 70
	progress = float(current)/float(total)

	if progress <= 1.0:
		sys.stdout.write('[')
		pos = int(barWidth*progress)
		for i in range(barWidth):
			if i < pos:
				sys.stdout.write('=')
			elif i == pos:
				sys.stdout.write('>')
			else:
				sys.stdout.write(' ')
		sys.stdout.write('] '+str(int(progress*100))+' %\r')
		sys.stdout.flush()

#--------------------------------------------------------------------------------------------------------------------------------------

#reference
reference = importReference(sys.argv[3])

#eventalign output
sixmer2eventsBrdU = {}
f = open(sys.argv[1],'r')
currentRead = ''
readCounter = 0
for line in f:

	splitLine = line.rstrip().split('\t')
	contig = splitLine[0]

	#ignore the header line
	if contig == 'contig':
		continue

	#ignore the complement for now
	if reverseComplement(splitLine[2]) == splitLine[9]:
		continue

	#ignore if the flanking bases are thymidines
	position = int(splitLine[1])
	if position < 1 or (position+5) >= len(reference[contig]):
		continue
	if reference[contig][position - 1] == 'T' or reference[contig][position - 2] == 'T' or reference[contig][position + 6] == 'T' or reference[contig][position + 5] == 'T':
		continue

	readIndex = splitLine[3]
	if readIndex != currentRead:
		currentRead = readIndex
		readCounter += 1
		displayProgress(readCounter, totalReads)

	if readCounter == totalReads:
		break

	sixmer = splitLine[9]
	if sixmer not in sixmer2eventsBrdU:
		sixmer2eventsBrdU[sixmer] = []
	else:
		sixmer2eventsBrdU[sixmer].append( float(splitLine[6]) )
f.close()

#pore model
model = {}
f = open(sys.argv[2],'r')
for line in f:

	if line[0] == '#' or line[0] == 'k':
		continue
	
	splitLine = line.rstrip().split('\t')
	model[splitLine[0]] = [	float(splitLine[1]), float(splitLine[2]) ]
f.close()


for i, key in enumerate(sixmer2eventsBrdU):

	if key == 'NNNNNN':
		continue
	
	#if key in mixture:
	if len( sixmer2eventsBrdU[key] ) > 40000:
		x = np.linspace( np.mean(sixmer2eventsBrdU[key])-15, np.mean(sixmer2eventsBrdU[key])+15, 1000 )
		plt.figure(i)

		#density
		densityBrdU = stats.kde.gaussian_kde( sixmer2eventsBrdU[key] )
		plt.plot(x, densityBrdU(x), label='BrdU Density')
		
		#pore model
		yModel = mlab.normpdf( x, model[key][0], model[key][1] )
		plt.plot(x, yModel, label='6mer Pore Model')

		#trimodal
		#yMix = mlab.normpdf( x, mixture[key][0], mixture[key][1] )
		#plt.plot( x, yMix, label='Fit Distribution (1)')
		#yMix = mlab.normpdf( x, mixture[key][2], mixture[key][3] )
		#plt.plot( x, yMix, label='Fit Distribution (2)')
		#yMix = mlab.normpdf( x, mixture[key][4], mixture[key][5] )
		#plt.plot( x, yMix, label='Fit Distribution (3)')

		#plotting stuff
		plt.xlabel('pA')
		plt.ylabel('Count')
		plt.title( key ) #+ '  N=' + str(int(mixture[key][6])) )
		plt.legend(loc='upper right')
		plt.savefig( key + '.png' )
		plt.close()
