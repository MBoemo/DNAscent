import matplotlib
matplotlib.use('Agg')
from matplotlib import pyplot as plt
import numpy as np
import sys

#Usage: python plotEventAlignments.py stderr.out DNAscent.detect
#where stderr.out is the result of a stderr redirect after running DNAscent detect with #define TEST_ALIGNMENT 1

maxReads = 2000
plot = False

threshold = 1.25 #log likelihood threshold for a positive BrdU call (this should be the default value used by DNAscent regions)

hist_emission = []
hist_gap = []
readID2emission = {}
readID2gap = {}
readID2residual = {}
readID2shift = {}
readID2scale = {}
readID2var = {}

x=[]
y=[]
progress = 0

#go through the alignments in order to:
# -plot the distribution of log emissions and gaps,
# -make individual plots of the alignments (if plot = True),
# -create a map from readID to emission so that we can crosscheck that against positive analogue call rates.
f = open(sys.argv[1],'r')
for line in f:

	if line[0] == '>':

		progress += 1
		if progress % 100 == 0:
			print(float(progress)/maxReads)
		if progress >= maxReads:
			break

		if len(x) > 0 and plot and progress < 10:
			plt.figure()
			plt.plot(x,y)
			plt.xlabel('Event Index')
			plt.ylabel('kmer Index')
			plt.savefig(readID+'.png')
			plt.close()

		if len(x) > 0:

			residual = np.polyfit(x,y,1,full=True)
			readID2residual[readID] = residual[1][0]/len(x)

		readID = line.rstrip()[1:]
		x=[]
		y=[]

	else:

		splitLine = line.rstrip().split()
		if splitLine[0] == 'avg_log_emission':
			hist_emission.append(float(splitLine[1]))
			readID2emission[readID] = float(splitLine[1])
		elif splitLine[0] == 'spanned':
			continue
		elif splitLine[0] == 'maxGap':
			hist_gap.append(int(splitLine[1]))
			readID2gap[readID] = float(splitLine[1])
		elif splitLine[0] == 'shift':
			readID2shift[readID] = float(splitLine[1])
		elif splitLine[0] == 'scale':
			readID2scale[readID] = float(splitLine[1])
		elif splitLine[0] == 'var':
			readID2var[readID] = float(splitLine[1])
		elif splitLine[0] == 'drift':
			continue
		else:	
			[eventIdx,kmerIdx] = splitLine
			x.append(int(eventIdx))
			y.append(int(kmerIdx))
f.close()

plt.figure()
plt.hist(hist_emission, 50)
plt.xlabel('Average Log Emission')
plt.ylabel('Count')
plt.savefig('hist_log_emission.pdf')
plt.close()

plt.figure()
plt.hist(hist_gap, 50)
plt.xlabel('Maximum Gap')
plt.ylabel('Count')
plt.savefig('hist_gap.pdf')
plt.close()

f = open(sys.argv[2],'r')
readID2callFraction = {}
first = True
progress = 0
for line in f:

	if line[0] == '>':

		progress += 1
		if progress % 100 == 0:
			print(float(progress)/maxReads)
		if progress >= maxReads:
			break

		if not first:
			readID2callFraction[readID] = float(brduCalls)/float(numAttempts)		

		else:
			first = False
		splitLine = line.rstrip().split()
		readID = splitLine[0][1:]

		numAttempts = 0
		brduCalls = 0
		methylCalls = 0
		brduMethyDeclined = 0

	else:

		splitLine = line.rstrip().split()
		logLikelihood = float(splitLine[1])
		if len(splitLine) > 4:
			logLikelihood_BrdUvsMethyl = float(splitLine[2])
			logLikelihood_MethylvsThym = float(splitLine[3])

			if logLikelihood > threshold and logLikelihood_BrdUvsMethyl < threshold:
				brduMethyDeclined += 1
			elif logLikelihood > threshold and logLikelihood_BrdUvsMethyl > threshold:
				brduCalls += 1
			elif logLikelihood_MethylvsThym > threshold and logLikelihood_BrdUvsMethyl < threshold:
				methylCalls += 1
			numAttempts += 1


		else:
			if logLikelihood > threshold:
				brduCalls += 1
			numAttempts += 1
f.close()

hist_x = []
hist_y = []

for ID in readID2callFraction:
	if ID in readID2emission:

		hist_x.append(readID2callFraction[ID])
		hist_y.append(readID2emission[ID])

plt.figure()
plt.scatter(hist_x,hist_y,alpha=0.5)
plt.xlabel('BrdU Calls (Calls/Attempts)')
plt.ylabel('Average Log Emission')
plt.savefig('scatter_emission_vs_calls.pdf')
plt.close()

hist_x = []
hist_y = []
for ID in readID2callFraction:
	if ID in readID2gap:

		hist_x.append(readID2callFraction[ID])
		hist_y.append(readID2gap[ID])

plt.figure()
plt.scatter(hist_x,hist_y,alpha=0.5)
plt.xlabel('BrdU Calls (Calls/Attempts)')
plt.ylabel('Max Gap')
plt.savefig('scatter_gap_vs_calls.pdf')
plt.close()

hist_x = []
hist_y = []
for ID in readID2callFraction:
	if ID in readID2residual:

		hist_x.append(readID2callFraction[ID])
		hist_y.append(readID2residual[ID])

plt.figure()
plt.scatter(hist_x,hist_y,alpha=0.5)
plt.xlabel('BrdU Calls (Calls/Attempts)')
plt.ylabel('Residual / Num Events')
plt.savefig('scatter_residual_vs_calls.pdf')
plt.close()

hist_x = []
hist_y = []
for ID in readID2callFraction:
	if ID in readID2shift:

		hist_x.append(readID2callFraction[ID])
		hist_y.append(readID2shift[ID])

plt.figure()
plt.scatter(hist_x,hist_y,alpha=0.5)
plt.xlabel('BrdU Calls (Calls/Attempts)')
plt.ylabel('Shift')
plt.savefig('scatter_shift_vs_calls.pdf')
plt.close()

hist_x = []
hist_y = []
for ID in readID2callFraction:
	if ID in readID2scale:

		hist_x.append(readID2callFraction[ID])
		hist_y.append(readID2scale[ID])

plt.figure()
plt.scatter(hist_x,hist_y,alpha=0.5)
plt.xlabel('BrdU Calls (Calls/Attempts)')
plt.ylabel('Scale')
plt.savefig('scatter_scale_vs_calls.pdf')
plt.close()

hist_x = []
hist_y = []
for ID in readID2callFraction:
	if ID in readID2var:

		hist_x.append(readID2callFraction[ID])
		hist_y.append(readID2var[ID])

plt.figure()
plt.scatter(hist_x,hist_y,alpha=0.5)
plt.xlabel('BrdU Calls (Calls/Attempts)')
plt.ylabel('Var')
plt.savefig('scatter_var_vs_calls.pdf')
plt.close()
