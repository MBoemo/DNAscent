#----------------------------------------------------------
# Copyright 2020 University of Cambridge
# Written by Michael A. Boemo (mb915@cam.ac.uk)
# This software is licensed under GPL-3.0.  You should have
# received a copy of the license with this software.  If
# not, please Email the author.
#----------------------------------------------------------

import sys
import os


#--------------------------------------------------------------------------------------------------------------------------------------
def splashHelp():
	s = """dnascent2bedgraph.py: Converts the output of DNAscent detect and forkSense into bedgraphs.
To run dnascent2bedgraph.py, do:
  python dnascent2bedgraph.py [arguments]
Example:
  python dnascent2bedgraph.py -d /path/to/dnascentDetect.out -f /path/to/dnascentForksense.out -o /path/to/newBedgraphDir
Required arguments are at least one of the following:
  -d,--detect               path to DNAscent detect output file,
  -f,--forkSense            path to DNAscent forkSense output file.
Required argument is:
  -o,--output               output directory which will be created.
Optional arguments are:
     --minLength            only convert reads with specified minimum read length (in base pairs) into bedgraphs (default: 1),
     --maxLength            only convert reads with specified maximum read length (in base pairs) into bedgraphs (default: Inf),
  -n,--maxReads             maximum number of reads to convert into bedgraphs (default: Inf),
     --targets              forkSense bed file with specific reads to plot,
     --filesPerDir          maximum reads per subdirectory (default: 300).
Written by Michael Boemo, Department of Pathology, University of Cambridge.
Please submit bug reports to GitHub Issues (https://github.com/MBoemo/DNAscent/issues)."""

	print(s)
	exit(0)


#--------------------------------------------------------------------------------------------------------------------------------------
class arguments:
	pass


#--------------------------------------------------------------------------------------------------------------------------------------
def parseArguments(args):

	a = arguments()

	a.minLength = 1
	a.maxLength = 1000000000
	a.maxReads = 1000000000
	a.filesPerDir = 300
	a.useTargets = False

	for i, argument in enumerate(args):
			
		if argument == '-d' or argument == '--detect':
			a.detectPath = str(args[i+1])

		elif argument == '-f' or argument == '--forkSense':
			a.sensePath = str(args[i+1])

		elif argument == '-o' or argument == '--output':
			a.outDir = str(args[i+1])

		elif argument == '--minLength':
			a.minLength = int(args[i+1])

		elif argument == '--maxLength':
			a.maxLength = int(args[i+1])

		elif argument == '-n' or argument == '--maxReads':
			a.maxReads = int(args[i+1])

		elif argument == '--targets':
			a.targetPath = str(args[i+1])
			a.useTargets = True

		elif argument == '--filesPerDir':
			a.filesPerDir = int(args[i+1])

		elif argument == '-h' or argument == '--help':
			splashHelp()

	#check that required arguments are met
	if not ( ( hasattr( a, 'detectPath') or hasattr( a, 'sensePath') ) and  hasattr( a, 'outDir') ):
		splashHelp() 
	return a


#--------------------------------------------------------------------------------------------------------------------------------------
def makeDetectLine(line, chromosome):
	splitLine = line.rstrip().split()
	pos = int(splitLine[0])
	probBrdU = float(splitLine[2])
	probEdU = float(splitLine[1])
	sixMer = splitLine[2]
	return (chromosome + ' ' + str(pos) + ' ' + str(pos+1) + ' ' + str(probBrdU) + '\n',chromosome + ' ' + str(pos) + ' ' + str(pos+1) + ' ' + str(probEdU) + '\n')


#--------------------------------------------------------------------------------------------------------------------------------------
def makeSenseLine(line, chromosome, prevPos):
	splitLine = line.rstrip().split()
	pos = int(splitLine[0])
	probEdUsegment = float(splitLine[1])
	probBrdUsegment = float(splitLine[2])
	return (chromosome + ' ' + str(prevPos) + ' ' + str(pos) + ' ' + str(probBrdUsegment) + '\n', chromosome + ' ' + str(prevPos) + ' ' + str(pos) + ' ' + str(probEdUsegment) + '\n', chromosome + ' ' + str(prevPos) + ' ' + str(pos) + ' ' + '\n')


#--------------------------------------------------------------------------------------------------------------------------------------
def parseBaseFile(fname, args, targetIDs):
	print('Parsing '+fname[0]+'...')
	first = True
	count = 0
	f = open(fname[0],'r')
	readID2directory = {}
	directoryCount = 0
	for line in f:

		#ignore blank lines just in case (but these shouldn't exist anyway)
		if not line.rstrip():
			continue

		#ignore the file header
		if line[0] == '#':
			continue
		
		if line[0] == '>':

			if not first:

				if args.useTargets and readID not in targetIDs:

					#get readID and chromosome
					splitLine = line.rstrip().split(' ')
					readID = splitLine[0][1:]
					chromosome = splitLine[1]
					strand = splitLine[4]
					mappingStart = int(splitLine[2])
					mappingEnd = int(splitLine[3])
					prevPos = mappingStart

					first = False
					buff = []

					continue

				rLen = mappingEnd - mappingStart

				if rLen > args.minLength and rLen < args.maxLength:

					if count % args.filesPerDir == 0:
						directoryCount += 1
						os.system('mkdir '+args.outDir + '/'+str(directoryCount))

					count += 1

					#stop if we've hit the max reads
					if count > args.maxReads:
						break

					readID2directory[readID] = directoryCount

					if fname[1] == "detect":
						f_bg = open( args.outDir + '/' + str(directoryCount) + '/' + readID + '.BrdUdetect.bedgraph','w')
						f_bg2 = open( args.outDir + '/' + str(directoryCount) + '/' + readID + '.EdUdetect.bedgraph','w')
						f_bg.write( 'track type=bedGraph name="'+readID +'" description="BedGraph format" visibility=full color=200,100,0 altColor=0,100,200 priority=20 viewLimits=0.0:1.0'+'\n')
						f_bg2.write( 'track type=bedGraph name="'+readID +'" description="BedGraph format" visibility=full color=93,197,186 altColor=0,100,200 priority=20 viewLimits=0.0:1.0'+'\n')

						for l in buff:
							f_bg.write(l[0])
							f_bg2.write(l[1])
						f_bg.close()
						f_bg2.close()

					elif fname[1] == "sense":
						#leftward moving fork
						f_BrdUsegment = open( args.outDir + '/' + str(readID2directory[readID]) + '/' + readID + '_BrdUsegment.bedgraph','w')
						f_BrdUsegment.write( 'track type=bedGraph name="'+readID + '_' + strand + '_BrdUsegment'+'" description="BedGraph format" visibility=full color=200,100,0 altColor=0,100,200 priority=20 viewLimits=0.0:1.0'+'\n')

						for l in buff:
							f_BrdUsegment.write(l[0])
						f_BrdUsegment.close()

						#rightward moving fork
						f_EdUsegment = open( args.outDir + '/' + str(readID2directory[readID]) + '/' + readID + '_EdUsegment.bedgraph','w')
						f_EdUsegment.write( 'track type=bedGraph name="'+readID + '_' + strand + '_EdUsegment'+'" description="BedGraph format" visibility=full color=93,197,186 altColor=0,100,200 priority=20 viewLimits=0.0:1.0'+'\n')

						for l in buff:
							f_EdUsegment.write(l[1])
						f_EdUsegment.close()
						
					
			#get readID and chromosome
			splitLine = line.rstrip().split(' ')
			readID = splitLine[0][1:]
			chromosome = splitLine[1]
			strand = splitLine[4]
			mappingStart = int(splitLine[2])
			mappingEnd = int(splitLine[3])
			prevPos = mappingStart

			first = False
			buff = []

		elif line[0] == '%':
			continue


		elif args.useTargets:

			if readID not in targetIDs:
				continue
			else:
				if fname[1] == "detect":
					buff.append( makeDetectLine(line,chromosome) )
				elif fname[1] == "sense":
					splitLine = line.rstrip().split()
					pos = int(splitLine[0])
					buff.append( makeSenseLine(line,chromosome,prevPos) )
					prevPos = pos
		else:

			if fname[1] == "detect":
				buff.append( makeDetectLine(line,chromosome) )
			elif fname[1] == "sense":
				splitLine = line.rstrip().split()
				pos = int(splitLine[0])
				buff.append( makeSenseLine(line,chromosome,prevPos) )
				prevPos = pos


	rLen = mappingEnd - mappingStart

	if rLen > args.minLength and rLen < args.maxLength and count < args.maxReads:

		if count % args.filesPerDir == 0:
			directoryCount += 1
			os.system('mkdir '+args.outDir + '/'+str(directoryCount))

		readID2directory[readID] = directoryCount

		if fname[1] == "detect":
			f_bg = open( args.outDir + '/' + str(directoryCount) + '/' + readID + '.BrdUdetect.bedgraph','w')
			f_bg2 = open( args.outDir + '/' + str(directoryCount) + '/' + readID + '.EdUdetect.bedgraph','w')
			f_bg.write( 'track type=bedGraph name="'+readID +'" description="BedGraph format" visibility=full color=200,100,0 altColor=0,100,200 priority=20 viewLimits=0.0:1.0'+'\n')
			f_bg2.write( 'track type=bedGraph name="'+readID +'" description="BedGraph format" visibility=full color=93,197,186 altColor=0,100,200 priority=20 viewLimits=0.0:1.0'+'\n')

			for l in buff:
				f_bg.write(l[0])
				f_bg2.write(l[1])
			f_bg.close()
			f_bg2.close()

		elif fname[1] == "sense":
			#leftward moving fork
			f_BrdUsegment = open( args.outDir + '/' + str(readID2directory[readID]) + '/' + readID + '_BrdUsegment.bedgraph','w')
			f_BrdUsegment.write( 'track type=bedGraph name="'+readID + '_' + strand + '_BrdUsegment'+'" description="BedGraph format" visibility=full color=200,100,0 altColor=0,100,200 priority=20 viewLimits=0.0:1.0'+'\n')

			for l in buff:
				f_BrdUsegment.write(l[0])
			f_BrdUsegment.close()

			#rightward moving fork
			f_EdUsegment = open( args.outDir + '/' + str(readID2directory[readID]) + '/' + readID + '_EdUsegment.bedgraph','w')
			f_EdUsegment.write( 'track type=bedGraph name="'+readID + '_' + strand + '_EdUsegment'+'" description="BedGraph format" visibility=full color=93,197,186 altColor=0,100,200 priority=20 viewLimits=0.0:1.0'+'\n')

			for l in buff:
				f_EdUsegment.write(l[1])
			f_EdUsegment.close()
			
	f.close()
	print('Done.')
	return readID2directory


#--------------------------------------------------------------------------------------------------------------------------------------
def parseSecondaryFile(fname, readID2directory,args, targetIDs):
	print('Parsing '+fname[0]+'...')
	f = open(fname[0],'r')
	first = True
	count = 0
	signalLineBuffer = []
	forkDirLineBuffer = []

	for line in f:

		if not line.rstrip():
			continue

		if line[0] == '#':
			continue
		
		if line[0] == '>':

			if not first:

				if args.useTargets and readID not in targetIDs:

					#get readID and chromosome
					splitLine = line.rstrip().split(' ')
					readID = splitLine[0][1:]
					chromosome = splitLine[1]
					strand = splitLine[-1:][0]
					mappingStart = int(splitLine[2])
					mappingEnd = int(splitLine[3])
					prevPos = mappingStart

					first = False
					buff = []

					continue

				rLen = mappingEnd - mappingStart

				if rLen > args.minLength and rLen < args.maxLength:

					count += 1

					if count > args.maxReads:
						break

					if readID in readID2directory:

						if fname[1] == "sense":
							#leftward moving fork
							f_BrdUsegment = open( args.outDir + '/' + str(readID2directory[readID]) + '/' + readID + '_BrdUsegment.bedgraph','w')
							f_BrdUsegment.write( 'track type=bedGraph name="'+readID + '_' + strand + '_BrdUsegment'+'" description="BedGraph format" visibility=full color=200,100,0 altColor=0,100,200 priority=20 viewLimits=0.0:1.0'+'\n')

							for l in buff:
								f_BrdUsegment.write(l[0])
							f_BrdUsegment.close()

							#rightward moving fork
							f_EdUsegment = open( args.outDir + '/' + str(readID2directory[readID]) + '/' + readID + '_EdUsegment.bedgraph','w')
							f_EdUsegment.write( 'track type=bedGraph name="'+readID + '_' + strand + '_EdUsegment'+'" description="BedGraph format" visibility=full color=93,197,186 altColor=0,100,200 priority=20 viewLimits=0.0:1.0'+'\n')

							for l in buff:
								f_EdUsegment.write(l[1])
							f_EdUsegment.close()
							
					
			#get readID and chromosome
			splitLine = line.rstrip().split(' ')
			readID = splitLine[0][1:]
			chromosome = splitLine[1]
			strand = splitLine[-1:][0]
			mappingStart = int(splitLine[2])
			mappingEnd = int(splitLine[3])
			prevPos = mappingStart

			first = False
			buff = []
		elif line[0] == '%':
			continue

		elif args.useTargets:

			if readID not in targetIDs:
				continue
			else:
				if fname[1] == "sense":
					splitLine = line.rstrip().split()
					pos = int(splitLine[0])
					buff.append( makeSenseLine(line,chromosome,prevPos) )
					prevPos = pos
		else:

			if fname[1] == "sense":
				splitLine = line.rstrip().split()
				pos = int(splitLine[0])
				buff.append( makeSenseLine(line,chromosome,prevPos) )
				prevPos = pos


	rLen = mappingEnd - mappingStart

	if rLen > args.minLength and rLen < args.maxLength and count < args.maxReads and readID in readID2directory:

		if fname[1] == "sense":
			#leftward moving fork
			f_BrdUsegment = open( args.outDir + '/' + str(readID2directory[readID]) + '/' + readID + '_BrdUsegment.bedgraph','w')
			f_BrdUsegment.write( 'track type=bedGraph name="'+readID + '_' + strand + '_BrdUsegment'+'" description="BedGraph format" visibility=full color=200,100,0 altColor=0,100,200 priority=20 viewLimits=0.0:1.0'+'\n')

			for l in buff:
				f_BrdUsegment.write(l[0])
			f_BrdUsegment.close()

			#rightward moving fork
			f_EdUsegment = open( args.outDir + '/' + str(readID2directory[readID]) + '/' + readID + '_EdUsegment.bedgraph','w')
			f_EdUsegment.write( 'track type=bedGraph name="'+readID + '_' + strand + '_EdUsegment'+'" description="BedGraph format" visibility=full color=93,197,186 altColor=0,100,200 priority=20 viewLimits=0.0:1.0'+'\n')

			for l in buff:
				f_EdUsegment.write(l[1])
			f_EdUsegment.close()
			

	print('Done.')
	f.close()


#--------------------------------------------------------------------------------------------------------------------------------------
#MAIN

args = parseArguments(sys.argv[1:])

targetIDs = []
if args.useTargets:
	f = open(args.targetPath,'r')
	for line in f:
		if line[0] == '#':
			continue
		splitLine = line.rstrip().split()
		readID = splitLine[3]
		targetIDs.append(readID)
	f.close()

#check the output 
if args.outDir[-1:] == "/":
	args.outDir = args.outDir[:-1]
if os.path.isdir(args.outDir):
	print('Output directory '+args.outDir+' already exists.  Exiting.')
	exit(0)
else:
	os.system('mkdir '+args.outDir)

baseFname = ""
secondaryFname = []

if hasattr( args, 'detectPath'):
	baseFname = (args.detectPath,"detect")
	if hasattr( args, 'sensePath'):
		secondaryFname.append((args.sensePath,"sense"))
else:
	if hasattr( args, 'sensePath'):
		baseFname = (args.sensePath,"sense")


readID2directory = parseBaseFile(baseFname, args, targetIDs)
for fname in secondaryFname:
	parseSecondaryFile(fname, readID2directory, args, targetIDs)
