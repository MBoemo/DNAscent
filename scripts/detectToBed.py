import sys

f = open(sys.argv[1],'r')
bed = open(sys.argv[1]+'.bed','w')

#write the header
bed.write('track name=detection description="Osiris-detected BrdU" useScore=1\n')

for line in f:
	if line[0] == '>':
		splitLine = line[1:].rstrip().split(' ')
		name = splitLine[0]
		splitLine2 = splitLine[1].split(':')
		chromosome = splitLine2[0]
		reset = True
		calls = 0
		attempts = 0

	else:
		splitLine = line.rstrip().split('\t')
		score = float(splitLine[1])
		pos = int(splitLine[0])
		
		attempts += 1
		if score > 1.5:
			calls += 1
		
		if reset:
			startPos = pos
			reset = False
		else:
			if pos - startPos >= 1000:
				reset = True
				bed.write(chromosome + ' ' + str(startPos) + ' ' + str(pos) + ' ' + chromosome + str(startPos) + ' ' + str(float(calls)/float(attempts)*1000*2) + '\n')
				attempts = 0
				calls = 0

bed.close()
f.close()
		
