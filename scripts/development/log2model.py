import sys
import numpy as np

f = open(sys.argv[1],'r')
f_out = open('all5mersFromLog.allModel','w')

dic = {}

for line in f:

	if line[0] == '>' or line[0] == 's':
		continue

	splitLine = line.rstrip().split('\t')
	if int(splitLine[0].split('_')[0]) == 16:
		splitStr = list(splitLine[1])
		splitStr[3] = 'B'
		splitStr = ''.join(splitStr)
		oriMu = float(splitLine[2])
		trMu = float(splitLine[3])
		oriSig = float(splitLine[4])
		trSig = float(splitLine[5])
		if splitStr in dic:
			dic[splitStr].append([trMu, trSig, oriMu, oriSig])
		else:
			dic[splitStr] = [[trMu, trSig, oriMu, oriSig]]

	elif int(splitLine[0].split('_')[0]) == 17:
		splitStr = list(splitLine[1])
		splitStr[2] = 'B'
		splitStr = ''.join(splitStr)
		oriMu = float(splitLine[2])
		trMu = float(splitLine[3])
		oriSig = float(splitLine[4])
		trSig = float(splitLine[5])
		if splitStr in dic:
			dic[splitStr].append([trMu, trSig, oriMu, oriSig])
		else:
			dic[splitStr] = [[trMu, trSig, oriMu, oriSig]]

	elif int(splitLine[0].split('_')[0]) == 18:
		splitStr = list(splitLine[1])
		splitStr[1] = 'B'
		splitStr = ''.join(splitStr)
		oriMu = float(splitLine[2])
		trMu = float(splitLine[3])
		oriSig = float(splitLine[4])
		trSig = float(splitLine[5])
		if splitStr in dic:
			dic[splitStr].append([trMu, trSig, oriMu, oriSig])
		else:
			dic[splitStr] = [[trMu, trSig, oriMu, oriSig]]

f_out.write('kmer\ttrMu\ttrSig\toriMu\toriSig\n')
for key in dic:
	means_bin = []
	for n in dic[key]:
		means_bin.append(n[0])
		f_out.write(key+'\t'+str(n[0])+'\t'+str(n[1])+'\t'+str(n[2])+'\t'+str(n[3])+'\n')
	f_out.write('>>\t'+str(np.mean(means_bin))+'\t'+str(np.std(means_bin))+'\n')
f.close()
f_out.close()
