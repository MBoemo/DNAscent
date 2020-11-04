import matplotlib
matplotlib.use('Agg')
from matplotlib import pyplot as plt
import sys
import numpy as np

#Usage: python callAgainstTcontent.py out.detect

maxReads = 1000

def reverseComplement(seq):
	newSeq = ''
	for s in seq:
		if s == 'A':
			newSeq += 'T'
		elif s == 'T':
			newSeq += 'A'
		elif s == 'C':
			newSeq += 'G'
		elif s == 'G':
			newSeq += 'C'
		else:
			warnings.warn("Nucleotides must be A, T, G, or C.")
			sys.exit()
			
	return newSeq[::-1]


#parse the detect file
thymCount2probs = {}
readCtr = 0
f = open(sys.argv[1],'r')
for line in f:
	if line[0] == '#':
		continue

	if line[0] == '>':
		[readID, chromosome, start, end, strand] = line.rstrip().split()

		readCtr += 1
		if readCtr > maxReads:
			break
	else:
		[pos, BrdUprob, kmerRef] = line.rstrip().split()

		if strand == 'rev':
			kmerRef = reverseComplement(kmerRef)

		key = kmerRef.count('T')

		if key not in thymCount2probs:
			thymCount2probs[key] = [float(BrdUprob)]
		else:
			thymCount2probs[key].append(float(BrdUprob))
f.close()

#reshape
x = []
y = []
y_err = []
for i in range(1,7):

	x.append(i)
	y.append(np.mean(thymCount2probs[i]))
	y_err.append(np.std(thymCount2probs[i]))

#plot
plt.figure()
plt.errorbar(x,y,yerr=y_err)
plt.xlabel('Number of Thymidines in Sixmer')
plt.ylabel('Average Called BrdU Probability')
plt.savefig('Calls_vs_Thymidines.pdf')
	
