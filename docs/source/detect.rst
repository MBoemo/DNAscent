.. _detect:

detect
===============================

``DNAscent detect`` is a ``DNAscent`` subprogram that analyses each nanopore-sequenced read and, at each thymidine position, assigns the probability that the base is really BrdU and EdU.

Usage
-----

.. code-block:: console

   To run DNAscent detect, do:
      DNAscent detect -b /path/to/alignment.bam -r /path/to/reference.fasta -i /path/to/index.dnascent -o /path/to/output.detect
   Required arguments are:
     -b,--bam                  path to alignment BAM file,
     -r,--reference            path to genome reference in fasta format,
     -i,--index                path to DNAscent index,
     -o,--output               path to output file that will be generated.
   Optional arguments are:
     -t,--threads              number of threads (default is 1 thread),
     --GPU                     use the GPU device indicated for prediction (default is CPU),
     -q,--quality              minimum mapping quality (default is 20),
     -l,--length               minimum read length in bp (default is 1000).

The main input of ``DNAscent detect`` is an alignment (bam file) between the sequence fastq from Guppy and the organism's reference genome.  This bam file should be sorted using ``samtools sort`` and indexed using ``samtools index`` so that there is a .bam.bai file in the same directory as the bam file. (Please see the example in :ref:`workflows` for details on how to do this.)  The full path to the reference genome used in the alignment should be passed using the ``-r`` flag, and the index required by the ``-i`` flag is the file created using ``DNAscent index`` (see :ref:`index_exe`).  

The number of threads is specified using the ``-t`` flag. ``DNAscent detect`` multithreads quite well by analysing a separate read on each thread, so multithreading is recommended. By default, the signal alignments and ResNet BrdU predictions are run on CPUs.  If a CUDA-compatible GPU device is specified using the ``--GPU`` flag, then the signal alignments will be run on CPUs using the threads specified with ``-t`` and the ResNet BrdU prediction will be run on the GPU. Your GPU device number can be found with the command ``nvidia-smi``. GPU use requires that CUDA and cuDNN are set up correctly on your system and that these libraries can be accessed. If they're not, DNAscent will default back to using CPUs.

It is sometimes useful to only run ``DNAscent detect`` on reads that exceed a certain mapping quality or length threshold (as measured by the subsequence of the contig that the read maps to).  In order to do this without having to filter the bam file, DNAscent provides the ``-l`` and ``-q`` flags.  Any read in the bam file with a reference length lower than the value specificed with ``-l`` or a mapping quality lower than the value specified with ``-q`` will be ignored.

Before calling BrdU and EdU in a read, ``DNAscent detect`` must first perform a fast event alignment (see https://www.biorxiv.org/content/10.1101/130633v2 for more details).  Quality control checks are performed on these alignments, and if they're not passed, then the read fails and is ignored.  Hence, the number of reads in the output file will be slightly lower than the number of input reads.  Typical failure rates are about 5-10%, although this will vary slightly depending on the read length, the BrdU substitution rate, and the genome sequenced.

Output
------

``DNAscent detect`` will produce a single human-readable output file with the name and location that you specified using the ``-o`` flag.  To aid organisation and reproducibility, each detect file starts with a short header.  The start of each header line is always a hash (#) character, and it specifies the input files and settings used, as well as the version and commit of DNAscent that produced the file.  An example is as follows:

.. code-block:: console

   #Alignment /path/to/alignment.bam
   #Genome /path/to/reference.fasta
   #Index /path/to/index.dnascent
   #Threads 1
   #Compute CPU
   #Mode CNN
   #MappingQuality 20
   #MappingLength 5000
   #SystemStartTime 09/02/2024 12:45:29
   #Software /path/to/DNAscent
   #Version 4.0.1
   #Commit 4cf80a7b89bdf510a91b54572f8f94d3daf9b167

You can easily access the header of any .detect file with ``head -11 /path/to/output.detect`` or, alternatively, ``grep '#' /path/to/output.detect``.

Below the header is data for each read.  Note that everything in this output file orients to the reference genome in the 5' --> 3' direction.  Each read starts with a line in the format:

.. code-block:: console

   >readID contig mappingStart mappingEnd strand

These lines always begin with a greater-than (>) character.  Therefore, an easy way to count the number of reads in the file is ``grep '>' detect.out | wc -l``.  The fields are:

* ``readID`` is a unique hexadecimal string assigned to each read by the Oxford Nanopore software,
* the read mapped between ``mappingStart`` and ``mappingEnd`` on ``contig``,
* ``strand`` either takes the value ``fwd``, indicating that the read mapped to the forward strand, or ``rev`` indicating that the read mapped to the reverse complement strand.

The following shows an example for a read that to the reverse strand between 48490 and 53033 on chromosome 1.

.. code-block:: console

   >0d64a203-81b5-4b6c-aa2f-67b20969a509 1 48490 53033 rev

Below these "start of read" lines, each line corresponds to the position of a thymidine in that read.  There are three tab-separated columns:

* the coordinate on the reference,
* probability that the thymidine is actually EdU,
* probability that the thymidine is actually BrdU


Consider the following examples:

.. code-block:: console

   >a4ea2872-9cb6-4218-afad-905f79204eb1 14 992440 996846 rev
   992448  0.125751        0.131483
   992450  0.082488        0.078428
   992451  0.070718        0.050604
   992453  0.062216        0.047409
   992456  0.056369        0.042582
   992457  0.046755        0.038603
   992459  0.056535        0.041545

   
Users may instead observe a "*" character in the fourth column; these characters mark indel sites for ``forkSense`` in order to avoid stalled fork false positives caused by indels in the genomic alignment. 


