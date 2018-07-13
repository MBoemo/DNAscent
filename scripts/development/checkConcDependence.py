import sys

percent_60 = {}
percent_80 = {}

f_80 = open(sys.argv[1],'r')
f_60 = open(sys.argv[2],'r')

for line in f_60:
	splitLine = line.rstrip().split('\t')
	percent_60[splitLine[0]] = splitLine[1:]
f_60.close()

for line in f_80:
	splitLine = line.rstrip().split('\t')
	sixMer = splitLine[0]

	BrdUguess_mean = float(splitLine[4])
	pi1 = float(splitLine[3])
	Thymguess_mean = float(splitLine[7])
	pi2 = float(splitLine[6])
	kl = float(splitLine[9])

	if kl > 2.5:
		if pi1 < float(percent_60[sixMer][2]):
			print sixMer, splitLine[1:]
			print sixMer, percent_60[sixMer]
	
