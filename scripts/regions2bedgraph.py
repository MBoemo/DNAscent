import sys
import os

printForkDir = False

f = open(sys.argv[1],'r')
first = True
printThisOne = False
count = 0
signalLineBuffer = []
forkDirLineBuffer = []

filesPerDir = 400
directoryCount = 0

for line in f:

	if not line.rstrip():
		continue
	
	if line[0] == '>':

		if not first:

			if printThisOne and len(signalLineBuffer) >= 5:

				if count % filesPerDir == 0:
					directoryCount += 1
					os.system('mkdir '+str(directoryCount))

				count += 1

				#signal
				outRegions = open( str(directoryCount) + '/' + str(count) + '_' + readID + '_scores.bedgraph','w')
				outRegions.write( 'track type=bedGraph name="'+readID + '_scores'+'" description="BedGraph format" visibility=full color=200,100,0 altColor=0,100,200 priority=20'+'\n')

				for l in signalLineBuffer:
					outRegions.write(l)
				outRegions.close()

				#fork direction
				if len( forkDirLineBuffer ) > 0 and printForkDir:

					outForkDir = open( str(directoryCount) + '/' + str(count) + '_' + readID + '_forkDir.bedgraph','w')
					outForkDir.write( 'track type=bedGraph name="'+readID + '_forkDir'+'" description="BedGraph format" visibility=full color=200,100,0 altColor=0,100,200 priority=20'+'\n')

					for l in forkDirLineBuffer:
						outForkDir.write(l)
					outForkDir.close()
				
		#get readID
		splitLine = line.rstrip().split(' ')
		readID = splitLine[0][1:]

		#get chromosome
		splitLine2 = splitLine[1].split(':')
		chromosome = splitLine2[0]

		first = False
		printThisOne = False
		signalLineBuffer = []
		forkDirLineBuffer = []

	else:
	
		splitLine = line.rstrip().split('\t')
		start = int(splitLine[0])
		end = int(splitLine[1])
		score = float(splitLine[2])
		call = splitLine[3]

		if call == "BrdU":
			printThisOne = True

		signalLineBuffer.append( chromosome + ' ' + str(start) + ' ' + str(end) + ' ' + str(score) + '\n' )
		

		if len(splitLine) > 4 and printForkDir:

			forkDir = splitLine[4]
			if forkDir == "-":

				forkDirLineBuffer.append( chromosome + ' ' + str(start) + ' ' + str(end) + ' ' + str(-5) + '\n' )
			else:

				forkDirLineBuffer.append( chromosome + ' ' + str(start) + ' ' + str(end) + ' ' + str(5) + '\n' )

f.close()

if printThisOne and len(signalLineBuffer) >= 5:
	count += 1
	#signal
	outRegions = open( str(directoryCount) + '/' +str(count) + '_' + readID + '_scores.bedgraph','w')
	outRegions.write( 'track type=bedGraph name="'+readID + '_scores'+'" description="BedGraph format" visibility=full color=200,100,0 altColor=0,100,200 priority=20'+'\n')

	for l in signalLineBuffer:
		outRegions.write(l)
	outRegions.close()

	#fork direction
	if len( forkDirLineBuffer ) > 0 and printForkDir:

		outForkDir = open( str(directoryCount) + '/' +str(count) + '_' + readID + '_forkDir.bedgraph','w')
		outForkDir.write( 'track type=bedGraph name="'+readID + '_forkDir'+'" description="BedGraph format" visibility=full color=200,100,0 altColor=0,100,200 priority=20'+'\n')

		for l in forkDirLineBuffer:
			outForkDir.write(l)
		outForkDir.close()
