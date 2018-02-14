#----------------------------------------------------------
# Copyright 2017 University of Oxford
# Written by Michael A. Boemo (michael.boemo@path.ox.ac.uk)
# This software is licensed under GPL-3.0.  You should have
# received a copy of the license with this software.  If
# not, please Email the author.
#----------------------------------------------------------

import pysam
import sys
import os
import gc
import h5py
import warnings
from joblib import Parallel, delayed
import multiprocessing

#--------------------------------------------------------------------------------------------------------------------------------------
class arguments:
	pass


#--------------------------------------------------------------------------------------------------------------------------------------
def splashHelp():
	s = """index.py: Osiris preprocessing script that will create an index of readIDs to fast5 file paths.
To run index.py, do:
  python index.py [arguments]
Example:
  python index.py -d /path/to/fast5 -t 20
Required arguments are:
  -d,--data                 path to top level directory of ONT fast5 reads.
Optional arguments are:
  -t,--threads              number of threads (default is 1 thread)."""

	print s
	exit(0)


#--------------------------------------------------------------------------------------------------------------------------------------
def parseArguments(args):

	a = arguments()
	a.threads = 1

	for i, argument in enumerate(args):
		if argument == '-t' or argument == '--threads':
			a.threads = int(args[i+1])

		elif argument == '-d' or argument == '--data':
			a.data = str(args[i+1])

		elif argument == '-h' or argument == '--help':
			splashHelp()
		elif argument[0] == '-':
			splashHelp()

	#check that required arguments are met
	if not hasattr( a, 'data'):
		splashHelp() 

	return a

#--------------------------------------------------------------------------------------------------------------------------------------
def get_readID(fast5file,root):
#opens a fast5 file, takes out the readID, and returns it

	if fast5file.endswith('.fast5'):
		
		try:
			#get signalID
			f_hdf5 = h5py.File(root+'/'+fast5file,'r')
			path = '/Raw/Reads'
			read_number = f_hdf5[path].keys()[0]
			readInfo = f_hdf5[path + '/' + read_number]
			readID = readInfo.attrs.get('read_id')
			return (readID, root+'/'+fast5file)

		except KeyError:
			warnings.warn('File '+root+'/'+fast5file+' did not have a valid readID path.  Skipping.', Warning)
			pass

		except IOError:
			warnings.warn('File '+root+'/'+fast5file+' could not be opened and may be corrupted.  Skipping.', Warning)


#MAIN--------------------------------------------------------------------------------------------------------------------------------------
args = sys.argv
a = parseArguments(args)

#output file to write on
fout = open('index.osiris','w')

progress = 0
#walk down the subdirectories and make a cumulative fasta in parallel
for root, dirs, files in os.walk(a.data, topdown=True):

	out = Parallel(n_jobs=a.threads)(delayed(get_readID)(f, root) for f in files)
	
	for r in out:
		readID, path = r
		fout.write(readID+'\t'+path+'\n')
	
	progress += len(out)
	sys.stdout.write("\rFinished indexing " + str(progress) + ' fast5 files...')
	sys.stdout.flush()

fout.close()
