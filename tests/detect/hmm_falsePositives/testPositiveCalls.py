import matplotlib
matplotlib.use('Agg')
from matplotlib import pyplot as plt
import numpy as np
import sys

#Usage: python testFalsePositives.py DNAscentDetect.stderr
#where DNAscentDetect.stderr is from stderr after running DNAscent detect with #define TEST_LL 1

llThreshold = 1.25

querySpan2positives = {}
querySpan2count = {}
events2positives = {}
events2count = {}

KL_scatter = []
ll_scatter = []

f = open(sys.argv[1],'r')
for ctr, line in enumerate(f):
	
	if line[0] == '<':

		if ctr > 1000000:
			break

		idx = 0
	else:
		idx += 1

		if idx == 1:
			KL = float(line.rstrip())
			KL_scatter.append(KL)
		if idx == 2:
			span = int(line.rstrip())				
		elif idx == 3:
			continue

		elif idx == 4:
			events = len(line.rstrip().split())
		elif idx == 5:
			ll = float(line.rstrip())
			ll_scatter.append(ll)
			#look at query
			if span in querySpan2positives:
				querySpan2count[span] = querySpan2count[span] + 1
				if ll > llThreshold:
					querySpan2positives[span] = querySpan2positives[span] + 1
			else:
				querySpan2count[span] = 1
				if ll > llThreshold:
					querySpan2positives[span] = 1
				else:
					querySpan2positives[span] = 0

			#look at event number
			if events in events2positives:
				events2count[events] += 1
				if ll > llThreshold:
					events2positives[events] += 1
			else:
				events2count[events] = 1
				if ll > llThreshold:
					events2positives[events] = 1
				else:
					events2positives[events] = 0

bar_x = []
bar_y = []
bar_z = []
#normalise
for query in querySpan2count:
	bar_x.append(query)
	bar_y.append(float(querySpan2positives[query])/querySpan2count[query])
	bar_z.append(querySpan2count[query])

plt.figure()
plt.bar(bar_x,bar_y)
plt.xlabel('Span on Query')
plt.ylabel('Positive Calls / Attempts')
plt.savefig('falsePositives_querySpan.pdf')
plt.close()

plt.figure()
plt.bar(bar_x,bar_z)
plt.xlabel('Span on Query')
plt.ylabel('Number of Call Positions')
plt.savefig('falsePositives_querySpan_positions.pdf')
plt.close()


bar_x = []
bar_y = []
bar_z = []
#normalise
for events in events2count:
	bar_x.append(events)
	bar_y.append(float(events2positives[events])/events2count[events])
	bar_z.append(events2count[events])

plt.figure()
plt.bar(bar_x,bar_y)
plt.xlabel('Number of Events')
plt.ylabel('Positive Calls / Attempts')
plt.savefig('falsePositives_eventNumbers.pdf')
plt.close()

plt.figure()
plt.bar(bar_x,bar_z)
plt.xlabel('Number of Events')
plt.ylabel('Number of Call Positions')
plt.savefig('falsePositives_eventNumbers_positions.pdf')
plt.close()

plt.figure()
plt.scatter(KL_scatter,ll_scatter,alpha=0.3)
plt.xlabel('Running KL Divergence')
plt.ylabel('Log Likelihood Ratio')
plt.savefig('falsePositives_runningKL.pdf')
plt.close()

