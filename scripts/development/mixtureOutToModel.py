import sys

f = open(sys.argv[1], 'r')
g = f.readlines()

print '#BrdU model trained by Osiris'
print 'kmer\tlevel_mean\tlevel_stdv'

for line in g[1:]:
	splitLine = line.split('\t')
	if float(splitLine[9]) > 2 and splitLine[0].find('T') != -1:
		print splitLine[0] + '\t' + splitLine[4] + '\t' + splitLine[5]
f.close()
