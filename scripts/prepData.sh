#!/bin/bash

#----------------------------------------------------------
# Copyright 2017 University of Oxford
# Written by Michael A. Boemo (michael.boemo@path.ox.ac.uk)
#----------------------------------------------------------

#Description:
#This shell script organises data from an Oxford Nanopore run for downstream analysis with Osiris.  It sorts reads by barcode, and performs a quality control on reads based on reference alignment.

#Dependencies:
# -pysam
# -samtools
# -bwa (will be installed if it doesn't exist)

#Troubleshooting:
# pysam > v0.9.1.2 (for version 0.7.7 there is an issue with the .reference_id attribute)
# There is a bug with Samtools version v1.19 that causes it to fall over when it's passed large files (it splits output into multiple files and then can't merge them back because it names them incorrectly).  This was fixed in > v1.2, but the Samtools sort -f flag has deprecated.  It's been replaced by -o.

##############
##  INPUTS  ##
##############

#directory for the FAST5 files produced in your run
FAST5_PASS_DIR=/HDD_archive/2016_11_30_hairpinBrdU_passOnly/pass

#reference file
REFERENCE=/data_disk_SSD/data/references/references_hairpin_withYcaps.fasta

#threads on which to run BWA alignment
THREADS=20


#################
##  DATA PREP  ##
#################

#check if bwa is compiled in the directory.  Clone the repository and compile it if it's not.
if [ ! -d bwa ]; then
	git clone https://github.com/lh3/bwa.git
	make -C bwa/
fi


#make symbolic links to the reference .fasta file, and create symbolic links to the FAST5 data files in directory data.fast5  
if [ ! -L reference.fasta ]; then 
	ln -s $REFERENCE reference.fasta
fi 

if [ ! -L data.fast5 ]; then
	ln -s $FAST5_PASS_DIR data.fast5
fi


#index the reference file with BWA
echo 'Indexing with BWA alignment...'
bwa index reference.fasta
echo 'Done.'

#create reads.fasta from data.fast5 using Osiris function import_2Dfasta


#do the alignment with BWA-MEM
echo 'Aligning...'
bwa mem -t $THREADS -k 1 -x ont2d reference.fasta reads.fasta | samtools view -Sb - | samtools sort -o alignments.sorted.bam - 
samtools index alignments.sorted.bam
echo 'Done.'


#produces one BAM file for each entry in the reference file
echo 'Assigning to construct...'
python assign_to_construct.py alignments.sorted.bam
echo 'Done.'


#index the construct bam files
echo 'Indexing BAM files...'
for bamfile in *.bam
do
	samtools index $bamfile
done
echo 'Done.'
