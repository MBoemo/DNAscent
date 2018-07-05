import sys
import math


#--------------------------------------------------------------------------------------------------------------------------------------
class arguments:
	pass


#--------------------------------------------------------------------------------------------------------------------------------------
def splashHelp():
	s = """callBrdU.py: Calls BrdU from a Nanopolish eventalign file.
To run callBrdU.py, do:
  python callBrdU.py [arguments]
Example:
  python callBrdU.py -tm /path/to/thymidine_model -bm /path/to/brdu_model -e /path/to/nanopolish_eventalignment -o output_prefix
Required arguments are:
  -tm,--ont_model           path to thymidine-only ONT model,
  -bm,--brdu_model          path to BrdU model,
  -e,--alignment            path to Nanopolish eventalign output file,
  -o,--output               output prefix.
Optional arguments are:
  -t,--threads              number of threads (default is 1 thread)."""

	print s
	exit(0)


#--------------------------------------------------------------------------------------------------------------------------------------
def importPoreModel( filename ):
	model = {}
	f = open(filename,'r')
	for line in f:

		if line[0] == '#' or line[0] == 'k':
			continue
	
		splitLine = line.rstrip().split('\t')
		model[splitLine[0]] = [	float(splitLine[1]), float(splitLine[2]) ]
	f.close()
	return model


#--------------------------------------------------------------------------------------------------------------------------------------
def evaluatePdf( mu, stdv, x ):
	
	return (1.0 / math.sqrt( 2.0*(math.pi)*math.pow(stdv,2.0)) ) * math.exp( - math.pow( x - mu, 2.0 ) / (2.0*math.pow(stdv, 2.0 ) ) )


#--------------------------------------------------------------------------------------------------------------------------------------
def parseArguments(args):

	a = arguments()
	a.threads = 1

	for i, argument in enumerate(args):
			
		if argument == '-tm' or argument == '--ont_model':
			a.ont_model = str(args[i+1])

		elif argument == '-bm' or argument == '--brdu_model':
			a.brdu_model = str(args[i+1])

		elif argument == '-e' or argument == '--alignment':
			a.eventalign = str(args[i+1])

		elif argument == '-o' or argument == '--output':
			a.outFile = str(args[i+1])

		elif argument == '-t' or argument == '--threads':
			a.threads = int(args[i+1])

		elif argument == '-h' or argument == '--help':
			splashHelp()

	#check that required arguments are met
	if not hasattr( a, 'ont_model') or not hasattr( a, 'brdu_model') or not hasattr( a, 'eventalign') or not hasattr( a, 'outFile'):
		splashHelp() 

	return a


#MAIN----------------------------------------------------------------------------------------------------------------------------------

#parse arguments
args = sys.argv
a = parseArguments(args)

#open output file
f_out = open(a.outFile + '.brdu','w')

#import models
brdu_model = importPoreModel( a.brdu_model )
ont_model = importPoreModel( a.ont_model )

#overall stats
BrdUPos = 0
ThymPos = 0
NaPos = 0

#iterate on eventalign lines
eventAlignment = open(a.eventalign,'r')
firstLine = True
currentRefPos = -1
llBuffer = []
linesBuffer = []
posChangeCounter = 0
numOfEventsCounter = 0
for line in eventAlignment:

	#skip the eventalign header line
	if firstLine:
		firstLine = False
		continue

	splitLine = line.split('\t')

	#get the position
	pos = int(splitLine[1])

	#if we changed positions, empty the buffer
	if pos != currentRefPos: 

		posChangeCounter += 1
		numOfEventsCounter += len(llBuffer)

		#average over buffer and give a hard call
		if len(llBuffer) > 0:
			avg = float(sum(llBuffer)) / float(len(llBuffer))
		else:
			avg = -1

		if avg >= 2.5:
			call = 1
			BrdUPos += 1
		elif avg == -1:
			call = 0
			NaPos += 1
		else:
			call = 0
			ThymPos += 1

		#write the lines
		for l in linesBuffer:
			f_out.write( l + str(avg) + '\t' + str(call) + '\n' )

		#empty the buffers 
		llBuffer = []
		linesBuffer = []

		#update current position
		currentRefPos = pos

		print BrdUPos, ThymPos, NaPos
		print float(numOfEventsCounter) / float(posChangeCounter)

	#compute log likelihood at this position
	sixMer = splitLine[9]
	eventMean = float( splitLine[6] )
	if sixMer.find('T') != -1 and sixMer in brdu_model:
		log_thym = math.log( evaluatePdf( ont_model[sixMer][0], ont_model[sixMer][1], eventMean ) )
		log_brdu = math.log( evaluatePdf( brdu_model[sixMer][0], brdu_model[sixMer][1], eventMean ) )
		llratio = log_brdu - log_thym
		llBuffer.append( llratio )
		linesBuffer.append( line.rstrip() + '\t' + str(llratio) + '\t' )
	else:
		linesBuffer.append( line.rstrip() + '\t' + str(-1) + '\t' )


eventAlignment.close()
f_out.close()




