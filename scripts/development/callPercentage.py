import sys


f = open(sys.argv[1],'r')
g = f.readlines()

calls = 0
tried = 0
for line in g:
	lineSplit = line.rstrip().split('\t')
	if float(lineSplit[1]) > 2.5:
		calls += 1

	if float(lineSplit[1]) != 0:
		tried += 1

print float(calls)/float(tried)
f.close()
