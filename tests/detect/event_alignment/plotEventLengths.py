import matplotlib
matplotlib.use('Agg')
from matplotlib import pyplot as plt
import sys

#Usage: python stderr.out
#where stderr.out is the result of a stderr redirect after running DNAscent detect #define EVENT_LENGTHS 1 in event_handling.cpp

lengths = []
f = open(sys.argv[1],'r')
for line in f:
	lengths.append(int(line.rstrip()))
f.close()

plt.figure()
plt.hist(lengths,25,log=True)
plt.xlabel('Event Length')
plt.ylabel('Count')
plt.savefig('eventLengths.pdf')
plt.close()
