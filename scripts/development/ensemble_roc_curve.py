import sys
import matplotlib
matplotlib.use('Agg')
from matplotlib import pyplot as plt

fn_pairs = [("/data/workspace/michael/2018_06_18_CAM_ONT_gDNA_BrdU_40_60_80_100/barcode11/detectionNew.threshold0.out","/data/workspace/michael/2018_06_18_CAM_ONT_gDNA_BrdU_40_60_80_100/barcode08/detectionNew.threshold0.out"), ("/data/workspace/michael/2018_06_18_CAM_ONT_gDNA_BrdU_40_60_80_100/barcode11/detectionNew.threshold1.0.out","/data/workspace/michael/2018_06_18_CAM_ONT_gDNA_BrdU_40_60_80_100/barcode08/detectionNew.threshold1.0.out"), ("/data/workspace/michael/2018_06_18_CAM_ONT_gDNA_BrdU_40_60_80_100/barcode11/detectionNew.threshold2.0.out","/data/workspace/michael/2018_06_18_CAM_ONT_gDNA_BrdU_40_60_80_100/barcode08/detectionNew.threshold2.0.out")]

plt.figure()
labels = ['6mer Inclusion Threshold = 0', '6mer Inclusion Threshold = 1.0', '6mer Inclusion Threshold = 1.5', '6mer Inclusion Threshold = 2.0']

for i, (fn_positive, fn_negative) in enumerate(fn_pairs):
	allThresholds = [x/float(10) for x in range(0,101,10)]
	xvalues = []
	yvalues = []

	for t in allThresholds:
		print fn_negative
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

		print fn_positive
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
	plt.plot(xvalues,yvalues,'o-',label=labels[i])
	for i, txt in enumerate(allThresholds[0:4]):
		plt.annotate(txt, (xvalues[i]+0.002,yvalues[i]-0.03))



plt.xlabel('False Positive Rate')
plt.ylabel('True Positive Rate')
plt.title('ROC Curve for Log Likelihood Thresholds')
plt.legend(loc='upper right', framealpha=0.5)
plt.savefig('ROC_curve_logLikelihood_ensemble.pdf')
