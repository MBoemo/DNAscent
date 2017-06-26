import pysam
import sys
import os
import gc
import h5py
import warnings

#--------------------------------------------------------------------------------------------------------------------------------------
class arguments:
	pass


#--------------------------------------------------------------------------------------------------------------------------------------
def splashHelp():
	s = """demultiplex.py: Osiris preprocessing script that will demultiplex barcoded reads.
To run demultiplex.py, do:
  python demultiplex.py [arguments]
Example:
  python demultiplex.py -r /path/to/reference.fasta -d /path/to/reads -t 20
Required arguments are:
  -r,--reference            path to reference file in fasta format,
  -d,--data                 path to top level directory of ONT reads,
Optional arguments are:
  -t,--threads              number of threads (default is 1 thread)."""

	print s
	exit(0)


#--------------------------------------------------------------------------------------------------------------------------------------
def parseArguments(args):

	a = arguments()

	for i, argument in enumerate(args):
		if argument == '-t' or argument == '--threads':
			a.threads = int(args[i+1])
			
		elif argument == '-r' or argument == '--reference':
			a.reference = str(args[i+1])

		elif argument == '-d' or argument == '--data':
			a.data = str(args[i+1])

		elif argument == '-h' or argument == '--help':
			splashHelp()

	return a


#--------------------------------------------------------------------------------------------------------------------------------------
def import_fasta(pathToReads, outFastaFilename):
#	takes a directory with fast5 nanopore reads at the top level, and extracts the 2D sequences in fasta format with the path to the file as the fasta header
#	ARGUMENTS
#       ---------
#	- pathToReads: full path to the directory that contains the fast5 files
#	  type: string
#	- outFastaFilename: filename for the output fasta file that contains all of the reads
#	  type: string
#	OUTPUTS
#       -------
#	- a fasta file written to the directory specified

	buffersize = 1024

	#output file to write on
	fout = open(outFastaFilename,'w')

	#path through the fast5 tree to get to the fastq sequence
	fast5path2fastq = '/Analyses/Basecall_1D_000/BaseCalled_template/Fastq'
	fast5path2score = '/Analyses/Basecall_1D_000/Summary/basecall_1d_template'

	#empty reads string, and count the number of subdirectories so we can print progress
	reads = ''
	numSubdirectories = len(next(os.walk(pathToReads, topdown=True))[1])
	readCount = 0

	#recursively go through the directory and subdirectories and extract fasta seqs until you reach the buffer, then write, release, and garbage collect
	for root, dirs, files in os.walk(pathToReads, topdown=True):

		for fast5file in files:

			readCount += 1

			if fast5file.endswith('.fast5'):
		
				#print progress every 5 subdirectories of reads
				if readCount % 10000 == 0:
					sys.stdout.write("\rExporting fast5 reads to fasta... read " + str(readCount))
					sys.stdout.flush()

				try:
					#open the fast5 file with h5py and grab the fastq
					ffast5 = h5py.File(root+'/'+fast5file,'r')
					fastq = ffast5[fast5path2fastq].value
					score = float(ffast5[fast5path2score].attrs.__getitem__('mean_qscore'))
					ffast5.close()
					fasta = fastq.split('\n')[1]
			
					#append the sequence in the fasta format, with the full path to the fast5 file as the sequence name
					if score > 10:
						reads += '>'+root+'/'+fast5file+'\n'+fasta+'\n'

				except KeyError:
					#warnings.warn('File '+root+'/'+fast5file+' did not have a valid fastq path.  Skipping.', Warning)
					pass

				except IOError:
					warnings.warn('File '+root+'/'+fast5file+' could not be opened and may be corrupted.  Skipping.', Warning)

				#write to the file and release the buffer
				if readCount % buffersize == 0:
					fout.write(reads)
					fout.flush()
					os.fsync(fout .fileno())
					reads = ''
					gc.collect()

		#flush the buffer and write once we're reached the end of fast5 files in the subdirectory
		fout.write(reads)
		fout.flush()
		os.fsync(fout .fileno())
		reads = ''
		gc.collect()
	
	#close output fasta file	
	fout.close()


#MAIN--------------------------------------------------------------------------------------------------------------------------------------
args = sys.argv
a = parseArguments(args)

import_fasta(a.data, os.getcwd()+'/reads.fasta')

os.system('bwa index ' + a.reference)
os.system('bwa mem -t '+str(a.threads)+' -x ont2d '+a.reference+' reads.fasta | samtools view -Sb - | samtools sort - alignments.sorted') 
os.system('samtools index alignments.sorted.bam')


sam_file = pysam.Samfile('alignments.sorted.bam')
out_files = list()

# open an output file for each reference sequence
for x in sam_file.references:
	print x
	out_files.append(pysam.Samfile(x + ".bam", "wb", template=sam_file))

for record in sam_file:
	ref_length = sam_file.lengths[record.reference_id]

	if record.aend is None or record.query_alignment_length is None:
		continue

	ref_cover = float(record.aend - record.pos) / ref_length
	query_cover = float(record.query_alignment_length) / record.query_length

	if ref_cover > 0.8 and query_cover > 0.8 and record.is_reverse == False:
		out_files[record.reference_id].write(record)
