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
def get_readID(fast5file,root):
#opens a fast5 file, takes out the readID, and returns it

	if fast5file.endswith('.fast5'):
		
		try:
			#get signalID
			fast5path2fastq = '/Analyses/Basecall_1D_000/BaseCalled_template/Fastq'
			ffast5 = h5py.File(root+'/'+fast5file,'r')
			fastq = ffast5[fast5path2fastq].value
			return fastq + '\n'

		except KeyError:
			warnings.warn('File '+root+'/'+fast5file+' did not have a valid readID path.  Skipping.', Warning)
			pass

		except IOError:
			warnings.warn('File '+root+'/'+fast5file+' could not be opened and may be corrupted.  Skipping.', Warning)


#MAIN--------------------------------------------------------------------------------------------------------------------------------------

#output file to write on
fout = open('reads.fastq','w')

progress = 0
#walk down the subdirectories and make a cumulative fasta in parallel
for root, dirs, files in os.walk(sys.argv[1], topdown=True):

	out = Parallel(n_jobs=50)(delayed(get_readID)(f, root) for f in files)
	
	for r in out:
		fout.write(r)
	
	progress += len(out)
	sys.stdout.write("\rFinished exporting " + str(progress) + ' fast5 files...')
	sys.stdout.flush()

fout.close()
