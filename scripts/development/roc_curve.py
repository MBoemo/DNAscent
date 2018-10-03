import sys
import matplotlib
matplotlib.use('Agg')
from matplotlib import pyplot as plt

fn_positive = "/data/workspace/michael/2018_06_18_CAM_ONT_gDNA_BrdU_40_60_80_100/barcode11/detectionTest.out"
fn_negative = "/data/workspace/michael/2018_06_18_CAM_ONT_gDNA_BrdU_40_60_80_100/barcode08/detectionTest.out"

allThresholds = [x/float(10) for x in range(0,101,10)]
xvalues = []
yvalues = []

for t in allThresholds:

	f = open(fn_negative,'r')
	numCalls = 0
	numAttempts = 0
	for line in f:

		if line[0] == '>':
			continue
		else:
			splitLine = line.split('\t')
			logLikelihood = float(splitLine[1])
			if logLikelihood > t:
				numCalls += 1
			numAttempts += 1
				
	falsePositives = float(numCalls)/float(numAttempts)
	xvalues.append(falsePositives)
	f.close()

	f = open(fn_positive,'r')
	numCalls = 0
	numAttempts = 0
	for line in f:

		if line[0] == '>':
			continue
		else:
			splitLine = line.split('\t')
			logLikelihood = float(splitLine[1])
			if logLikelihood > t:
				numCalls += 1
			numAttempts += 1
				
	truePositives = 2*(float(numCalls)/float(numAttempts)) #multiplied by two because maximum is 50%
	yvalues.append(truePositives)

	f.close()
plt.figure()
plt.plot(xvalues,yvalues,'o-')

for i, txt in enumerate(allThresholds[0:11]):
	plt.annotate(txt, (xvalues[i],yvalues[i]), (xvalues[i]+0.002,yvalues[i]-0.03) )

plt.xlabel('False Positive Rate')
plt.ylabel('True Positive Rate')
plt.xlim(0,.1)
plt.ylim(0,1)
plt.title('ROC Curve for Log Likelihood Threshold')
plt.savefig('ROC_curve_logLikelihood.pdf')
