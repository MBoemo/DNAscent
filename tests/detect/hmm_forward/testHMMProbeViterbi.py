from pomegranate import *
import matplotlib
matplotlib.use('Agg')
from matplotlib import pyplot as plt

#Usage: python testHMMForward.py DNAscentDetect.stderr
#where DNAscentDetect.stderr is from stderr after running DNAscent detect with #define TEST_HMM 1

def import_poreModel(filename):
#	takes the filename of an ONT pore model file and returns a map from kmer (string) to [mean,std] (list of floats)
#	ARGUMENTS
#       ---------
#	- filename: path to an ONT model file
#	  type: string
#	OUTPUTS
#       -------
#	- kmer2MeanStd: a map, keyed by a kmer, that returns the model mean and standard deviation signal for that kmer
#	  type: dictionary

	f = open(filename,'r')
	g = f.readlines()
	f.close()

	kmer2MeanStd = {}
	for line in g:
		if line[0] != '#' and line[0:4] != 'kmer': #ignore the header
			splitLine = line.split('\t')
			kmer2MeanStd[ splitLine[0] ] = [ float(splitLine[1]), float(splitLine[2]) ]
	g = None

	return kmer2MeanStd


def build_TrainingHMM(refSequence, thymidineModel, brduModel, scalings, useBrdU):

	hmm = HiddenMarkovModel()

	windowSize = 10
	brduStart = windowSize - 5
	brduEnd = windowSize + refSequence[windowSize:windowSize+6].rfind('T')

	refLength = len(refSequence)

	#new HMM transition parameters
	internalM12I = 0.3475
	internalI2I = 0.5
	internalM12M1 = 0.4

	externalD2D = 0.3
	externalD2M1 = 0.7
	externalI2M1 = 0.5
	externalM12D = 0.0025
	externalM12M1 = 0.25

	###########################################################################################################################
	# Add States to Model
	#Create the HMM states.  Iterate through the reference sequence, and make the repeating HMM module for each position in the sequence.

	#two-dimensional array for the states, where the columns are the positions in the reference
	states = [[0 for x in range(refLength)] for y in range(3)]

	#the first base
	i = 0
	emissions = thymidineModel[refSequence[0:6]]
	level_mu = scalings[0] + scalings[1] * emissions[0]
	level_sig = scalings[2] * emissions[1]

	#print(refSequence[i:i+6])
	#print(level_mu,level_sig)

	states[0][i] = State( UniformDistribution(0, 250, frozen=True), name='Insertion_'+str(i) )    
	states[1][i] = State( NormalDistribution(level_mu, level_sig), name='Match_'+str(i) )    
	states[2][i] = State( None, name='Deletion_'+str(i) )  
	for j in range(3):
		hmm.add_state(states[j][i])

	#make an insertion state before the first base
	firstI = State( UniformDistribution(0, 250, frozen=True), name='Insertion_Pre' )
	hmm.add_state(firstI)

	#insertion state before the first base
	hmm.add_transition(hmm.start, firstI, 0.25) #start to the first insertion
	hmm.add_transition(firstI, firstI, 0.25) #self loop

	#to the base 1 insertion
	hmm.add_transition(states[0][i], states[0][i], internalI2I , group='internal_I-to-I')
	hmm.add_transition(states[1][i], states[0][i], internalM12I , group='internal_M-to-I')

	#to the base 1 match
	hmm.add_transition(firstI, states[1][i], 0.5) #first insertion to first match
	hmm.add_transition(states[1][i], states[1][i], internalM12M1 , group='internal_M-to-M')
	hmm.add_transition(hmm.start, states[1][0], 0.5)  #start to M
	
	#to the base 1 deletion
	hmm.add_transition(firstI, states[2][i], 0.25) #first insertion to the first deletion
	hmm.add_transition(hmm.start, states[2][0], 0.25) #start to D

	#the rest of the sequence
	for i, char in enumerate(refSequence[1:-5]):

		i += 1

		if (useBrdU and brduStart <= i and i <= brduEnd):
			emissions = brduModel[refSequence[i:i+6]]
		else:
			emissions = thymidineModel[refSequence[i:i+6]]

		#correct for shift/scale/var
		level_mu = scalings[0] + scalings[1] * emissions[0]
		level_sig = scalings[2] * emissions[1]

		#print(refSequence[i:i+6])
		#print(level_mu,level_sig)
		
		#create states for this nucleotide
		states[0][i] = State( UniformDistribution(0, 250, frozen=True), name='Insertion_'+str(i) )    
		states[1][i] = State( NormalDistribution(level_mu, level_sig), name='Match_'+str(i) )    
		states[2][i] = State( None, name='Deletion_'+str(i) )     

		#add the state objects to the hmm model object
		for j in range(3):
			hmm.add_state(states[j][i])

		#internal transitions for this nucleotide
		hmm.add_transition(states[1][i], states[1][i], internalM12M1 , group='internal_M-to-M')
		hmm.add_transition(states[1][i], states[0][i], internalM12I , group='internal_M-to-I')
		hmm.add_transition(states[0][i], states[0][i], internalI2I , group='internal_I-to-I')

	#this is really just a safety thing... get the last index in the iterator
	for last_i,x in enumerate(refSequence[:-5]): 
		pass	

	###########################################################################################################################
	# Handle Transitions Between Modules, Handle Analogue Branch

	#We have to reach forward to the next state (state i+1) so it's easier to just do this in a separate loop from the internal transitions one
	for i, char in enumerate(refSequence[:-5]): 

		#Don't execute this if we're at the end, because there's no i+1 to reach forward to.
		if i != last_i:

			hmm.add_transition(states[1][i], states[2][i+1], externalM12D, group='external_M-to-D')
			hmm.add_transition(states[2][i], states[2][i+1], externalD2D, group='external_D-to-D')
			hmm.add_transition(states[2][i], states[1][i+1], externalD2M1, group='external_D-to-M')
			hmm.add_transition(states[1][i], states[1][i+1], externalM12M1, group='external_M-to-M')
			hmm.add_transition(states[0][i], states[1][i+1], externalI2M1, group='external_I-to-M')

	###########################################################################################################################
	# Handle Start and End

	#handle end states
	hmm.add_transition(states[0][last_i], hmm.end, externalI2M1 )
	hmm.add_transition(states[1][last_i], hmm.end, externalM12M1 + externalM12D)
	hmm.add_transition(states[2][last_i], hmm.end, 1.0)

	#bake the model
	#hmm.bake(merge='all',verbose=True)
	hmm.bake()
	return hmm

###########################################################################################################################
# MAIN

verbose = False

print('Loading pore models...')
f_thymidineModel = '/home/mb915/rds/hpc-work/development/DNAscent_dev/pore_models/template_median68pA.6mer.model'
thymidineModel = import_poreModel(f_thymidineModel)
f_brduModel = '/home/mb915/rds/hpc-work/development/DNAscent_dev/pore_models/BrdU.model'
brduModel = import_poreModel(f_brduModel)
print('Done.')

#scalings: shift,scale,var

allRatios = []
brduProb = []
thymProb = []

brduDel = []
brduIns = []
thymDel = []
thymIns = []

maxEvents = 10000

f = open(sys.argv[1],'r')
idx = 0
eventsRead = 0
for line in f:
	idx += 1
	if line[0] == '<':
		idx = 0
	if idx == 1:
		useBrdU = int(line.rstrip())
	elif idx == 2:
		scalings = line.rstrip().split()
		scalings = [float(i) for i in scalings]	
	elif idx == 3:
		sequence = line.rstrip()
	elif idx == 4:
		events = line.rstrip().split()
		events = [float(i) for i in events]
	elif idx == 5:
		eventsRead += 1
		hmm = build_TrainingHMM(sequence, thymidineModel, brduModel, scalings, useBrdU)
		vitOut = hmm.viterbi(events)
		viterbiStates = vitOut[1]
		names = []
		insertions = 0
		deletions = 0
		for s in viterbiStates:
			names.append(s[1].name)
			if 'Insertion' in s[1].name:
				insertions += 1
			elif 'Deletion' in s[1].name:
				deletions += 1
		if useBrdU:
			DNAscentProb_BrdU = float(line.rstrip())
			pomegranateProb_BrdU = hmm.log_probability(events)
			namesBrdU = names
			insertions_BrdU = insertions
			deletions_BrdU = deletions
		else:
			DNAscentProb_thym = float(line.rstrip())
			pomegranateProb_thym = hmm.log_probability(events)
			namesThym = names
			insertions_Thym = insertions
			deletions_Thym = deletions
			if verbose:
				print('------------------------------------------')
				print('Log Likelihood Ratio: ',DNAscentProb_BrdU-DNAscentProb_thym)
				print('Thymidine Insertions:',insertions_Thym)
				print('Thymidine Deletions:',deletions_Thym)
				print('BrdU Insertions:',insertions_BrdU)
				print('BrdU Deletions:',deletions_BrdU)
				print('Thym Path:',namesThym)
				print('BrdU Path:',namesBrdU)
				print('------------------------------------------')

			allRatios.append(DNAscentProb_BrdU-DNAscentProb_thym)
			brduProb.append(DNAscentProb_BrdU)
			thymProb.append(DNAscentProb_thym)

			brduIns.append(insertions_BrdU)
			brduDel.append(deletions_BrdU)
			thymIns.append(insertions_Thym)
			thymDel.append(deletions_Thym)


			if eventsRead % 1000:
				print(eventsRead/float(maxEvents))

			if eventsRead > maxEvents:
				break

f.close()


plt.figure()
plt.scatter(allRatios,brduProb,alpha=0.3)
plt.scatter(allRatios,thymProb,alpha=0.3)
plt.legend(['BrdU HMM','Thymidine HMM'])
plt.ylabel('Log Likelihood')
plt.xlabel('Log Likelihood Ratio (BrdU-to-Thym)')
plt.savefig('test_llRatios_Probability.pdf')
plt.close()

plt.figure()
plt.scatter(allRatios,brduIns,alpha=0.3)
plt.scatter(allRatios,thymIns,alpha=0.3)
plt.legend(['BrdU HMM','Thymidine HMM'])
plt.ylabel('Number of Viterbi Insertions')
plt.xlabel('Log Likelihood Ratio (BrdU-to-Thym)')
plt.savefig('test_llRatios_Insertions.pdf')
plt.close()

plt.figure()
plt.scatter(allRatios,brduDel,alpha=0.3)
plt.scatter(allRatios,thymDel,alpha=0.3)
plt.legend(['BrdU HMM','Thymidine HMM'])
plt.ylabel('Number of Viterbi Deletions')
plt.xlabel('Log Likelihood Ratio (BrdU-to-Thym)')
plt.savefig('test_llRatios_Deletions.pdf')
plt.close()


