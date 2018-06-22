import sys

#Useage: python findGenomicContexts.py [reference] [6mer]
#Example: python findGenomicContexts.py reference.fasta AGCCGA

#--------------------------------------------------------------------------------------------------------------------------------------
def reverseComplement(sequence):

	#seq to upper case
	sequence = sequence.upper()

	#build the complement
	revComp = ''
	for char in sequence:
		if char == 'A':
			revComp += 'T'
		elif char == 'T':
			revComp += 'A'
		elif char == 'C':
			revComp += 'G'
		elif char == 'G':
			revComp += 'C'
		#guard against illegal characters
		else:
			warnings.warn('Warning: Illegal character in sequence.  Legal characters are A, T, G, and C.', Warning)

	#take the reverse of the complement and return it
	return revComp[::-1]


#--------------------------------------------------------------------------------------------------------------------------------------
def importReference(filename):

	f = open(filename,'r')
	g = f.readlines()
	f.close()

	reference = {}
	for line in g:

		if line[0] == '>':
			currentChromosome = line[1:].rstrip()
			reference[ currentChromosome ] = ''

		if line[0] != '>':
			reference[ currentChromosome ] += line.rstrip().upper()

	return reference

#MAIN----------------------------------------------------------------------------------------------------------------------------------
reference = importReference( sys.argv[1] )
target6mer = sys.argv[2]

threePrimeFlanking = {'A':0, 'T':0, 'G':0, 'C':0}
fivePrimeFlanking = {'A':0, 'T':0, 'G':0, 'C':0}
total = 0

#template
for key in reference:

	if key == 'chrM':
		continue

	loc = 0
	start = 0
	while loc != -1:
		loc = reference[key].find(target6mer, start, len(reference[key]))
		start = loc + 1

		if loc != -1 and not (key == 'chrXII' and loc > 450000 and loc < 470000):
			print 'template', key, loc, reference[key][max(0,loc-5):min(len(reference[key]), loc+11)]
			total += 1
			if loc-1 > 0:
				fivePrimeFlanking[ reference[key][loc-1] ] += 1
			if loc+6 < len(reference[key]):
				threePrimeFlanking[ reference[key][loc+6] ] += 1

#complement
for key in reference:

	if key == 'chrM':
		continue

	refRC = reverseComplement(reference[key])

	loc = 0
	start = 0
	while loc != -1:
		loc = refRC.find(target6mer, start, len(refRC))
		start = loc + 1

		if loc != -1 and not (key == 'chrXII' and loc > (len(refRC) - 470000) and loc < (len(refRC) - 450000)):
			print 'complement', key, loc, refRC[max(0,loc-5):min(len(refRC), loc+11)]
			total += 1
			if loc-1 > 0:
				fivePrimeFlanking[ refRC[loc-1] ] += 1
			if loc+6 < len(refRC):
				threePrimeFlanking[ refRC[loc+6] ] += 1

#print statistics
print 'Total hits:', total
print '5prime Flanking Base (Fraction):'
for key in fivePrimeFlanking:
	print key, float(fivePrimeFlanking[key])/float(total)
print '3prime Flanking Base (Fraction):'
for key in threePrimeFlanking:
	print key, float(threePrimeFlanking[key])/float(total)

		
