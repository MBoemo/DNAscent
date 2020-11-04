import matplotlib
matplotlib.use('Agg')
from matplotlib import pyplot as plt
import sys

#Usage: python callsAgainstKL.py out.detect DNAscent/pore_models/BrdU_full_noThreshold.model


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


callThreshold = 2.5
#parse the detect file
kmer2calls = {}
kmer2attempts = {}
f = open(sys.argv[1],'r')
for line in f:
	if line[0] == '>':
		[readID, chromosome, start, end, strand] = line.rstrip().split()
	else:
		[pos, score, kmerRef, kmerRead] = line.rstrip().split()

		if strand == 'rev':
			kmerRef = reverseComplement(kmerRef)

		if kmerRef not in kmer2attempts:
			kmer2attempts[kmerRef] = 1
		else:
			kmer2attempts[kmerRef] += 1

		score = float(score)
		if score < callThreshold:
			continue

		if kmerRef not in kmer2calls:
			kmer2calls[kmerRef] = 1
		else:
			kmer2calls[kmerRef] += 1
f.close()

#parse the model file
kmer2KL = {}
f = open(sys.argv[2],'r')
for line in f:
	splitLine = line.rstrip().split()
	KL = float(splitLine[-1:][0])
	kmer = splitLine[0]
	kmer2KL[kmer] = KL
f.close()

#reshape
x = []
y = []
for kmer in kmer2calls:
	if kmer in kmer2KL:
		x.append(float(kmer2calls[kmer])/kmer2attempts[kmer])
		y.append(kmer2KL[kmer])

#plot
plt.figure()
plt.scatter(x,y,alpha=0.2)
plt.xlabel('(Positive Calls)/(Attempts)')
plt.ylabel('KL-Divergence')
plt.savefig('Calls_vs_Divergence.pdf')
	
