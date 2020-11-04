.. _cookbook:

Python Cookbook
===============================

The output file formats of all DNAscent executables were specifically designed to be easy to parse with short (Python, Perl, etc.) scripts with the aim of making it simple for users to make application-specific plots.  Here, we provide a brief "cookbook" of barebones analysis scripts that can be copied and modified by users.

The following barebones script parses the output of ``DNAscent detect``.  We iterate line-by-line and parse each field in the file.  

.. code-block:: python

   f = open('path/to/output.detect','r')

   for line in f:

	#ignore the header lines
   	if line[0] == '#':
		continue
	
	#split the line into a list by whitespace
	splitLine = line.rstrip().split()

	if line[0] == '>':

		readID = splitLine[0][1:]
		chromosome = splitLine[1]
		refStart = int(splitLine[2])
		refEnd = int(splitLine[3])
		strand = splitLine[4]
	else:
		posOnRef = int(splitLine[0])
		probBrdU = float(splitLine[1])
		sixMerOnRef = splitLine[2]

		#add these values to a container or do some processing here

   f.close()

The following barebones script parses the output of ``DNAscent forkSense``.  Note the similarity to the above script: All DNAscent output files were designed to have a very similar format to aid user processing.  

.. code-block:: python

   f = open('path/to/output.forkSense','r')

   for line in f:

	#ignore the header lines
   	if line[0] == '#':
		continue
	
	#split the line into a list by whitespace
	splitLine = line.rstrip().split()

	if line[0] == '>':

		readID = splitLine[0][1:]
		chromosome = splitLine[1]
		refStart = int(splitLine[2])
		refEnd = int(splitLine[3])
		strand = splitLine[4]
	else:
		posOnRef = int(splitLine[0])
		probLeftFork = float(splitLine[1])
		probRightFork = float(splitLine[2])

		#add these values to a container or do some processing here

   f.close()

And again for ``DNAscent regions``:

.. code-block:: python

   f = open('path/to/output.regions','r')

   for line in f:

	#ignore the header lines
   	if line[0] == '#':
		continue
	
	#split the line into a list by whitespace
	splitLine = line.rstrip().split()

	if line[0] == '>':

		readID = splitLine[0][1:]
		chromosome = splitLine[1]
		refStart = int(splitLine[2])
		refEnd = int(splitLine[3])
		strand = splitLine[4]
	else:
		regionStart = int(splitLine[0])
		regionEnd = int(splitLine[1])
		regionScore = float(splitLine[2])

		#add these values to a container or do some processing here

   f.close()


