import sys

f = open(sys.argv[1],'r')
g = f.readlines()

toSort = []
for line in g[1:]:
	toSort.append( line )

f.close()

toSort = sorted(toSort)

f = open(sys.argv[1]+'.sorted','w')
f.write(g[0])
for s in toSort:
	f.write(s)
f.close()
