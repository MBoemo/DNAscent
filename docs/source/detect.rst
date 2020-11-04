.. _detect:

Detect
===============================

``DNAscent detect`` is a ``DNAscent`` subprogram that goes through each read and, at each thymidine position, assigns the probability the thymidine is BrdU.

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
     --GPU                      use the GPU device indicated for prediction (default is CPU),
     -q,--quality              minimum mapping quality (default is 20),
     -l,--length               minimum read length in bp (default is 100).

The main input of ``DNAscent detect`` is an alignment (bam file) between the sequence fastq from Guppy and the organism's reference genome.  This bam file should be sorted using ``samtools sort`` and indexed using ``samtools index`` so that there is a .bam.bai file in the same directory as the bam file. (Please see the example in :ref:`workflows` for details on how to do this.)  The full path to the reference genome used in the alignment should be passed using the ``-r`` flag, and the index required by the ``-i`` flag is the file created using ``DNAscent index`` (see :ref:`index`).  

The number of threads is specified using the ``-t`` flag. ``DNAscent detect`` multithreads quite well by analysing a separate read on each thread, so multithreading is recommended. By default, the signal alignments and ResNet BrdU predictions are run on CPUs.  If a CUDA-compatible GPU device is specified using the ``--GPU`` flag, then the signal alignments will be run on CPUs using the threads specified with ``-t`` and the ResNet BrdU prediction will be run on the GPU. Your GPU device number can be found with the command ``nvidia-smi``. GPU use requires that CUDA and cuDNN are set up correctly on your system and that these libraries can be accessed. If they're not, DNAscent will default back to using CPUs.

It is sometimes useful to only run ``DNAscent detect`` on reads that exceed a certain mapping quality or length threshold (as measured by the subsequence of the contig that the read maps to).  In order to do this without having to filter the bam file, DNAscent provides the ``-l`` and ``-q`` flags.  Any read in the bam file with a reference length lower than the value specificed with ``-l`` or a mapping quality lower than the value specified with ``-q`` will be ignored.

Before calling BrdU in a read, ``DNAscent detect`` must first perform a fast event alignment (see https://www.biorxiv.org/content/10.1101/130633v2 for more details).  Quality control checks are performed on these alignments, and if they're not passed, then the read fails and is ignored.  Hence, the number of reads in the output file will be slightly lower than the number of input reads.  Typical failure rates are about 5-10%, although this will vary slightly depending on the read length, the BrdU substitution rate, and the genome sequenced.

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
   #SignalDilation 1.000000
   #Version 2.0.0
   #Commit 4cf80a7b89bdf510a91b54572f8f94d3daf9b167

You can easily access the header of any .detect file with ``head -11 /path/to/output.detect`` or, alternatively, ``grep '#' /path/to/output.detect``.

Below the header is data for each read.  Note that everything in this output file orients to the reference genome in the 5' --> 3' direction.  Each read starts with a line in the format:

.. code-block:: console

   >readID contig mappingStart mappingEnd strand

These lines always begin with a greater-than (>) character.  Therefore, an easy way to count the number of reads in the file is ``grep '>' detect.out | wc -l``.  The fields are:

* ``readID`` is a unique hexadecimal string assigned to each read by the Oxford Nanopore software,
* the read mapped between ``mappingStart`` and ``mappingEnd`` on ``contig``,
* ``strand`` either takes the value ``fwd``, indicating that the read mapped to the forward strand, or ``rev`` indicating that the read mapped to the reverse complement strand.

The following shows an example for a read that to the reverse strand between 239248 and 286543 on chrII.

.. code-block:: console

   >c602f23f-e892-42ba-8140-da949abafbdd chrII 239248 286543 rev

Below these "start of read" lines, each line corresponds to the position of a thymidine in that read.  There are three tab-separated columns:

* the coordinate on the reference,
* probability that the thymidine is actually BrdU,
* 6mer on the reference.


Consider the following examples:

.. code-block:: console

   >c6785e1f-10d2-49cb-8ca3-e8d48979001b chrXIII 74003 81176 rev
   74010   0.012874        TCTCTA
   74011   0.012428        CTCTAA
   74014   0.016811        TAACGA
   74017   0.013372        CGACCA
   74018   0.013836        GACCAA

Here, we're looking at the sequence TCTCTAACGACCAA on the reference genome.  Because this read maps to the reverse complement, a call is made at every A (instead of T) on the reference.  If instead we looked at a read that mapped to the forward strand, an example would be:

.. code-block:: console

   >5d10eb9a-aae1-4db8-8ec6-7ebb34d32575 chrXIII 72319 77137 fwd
   72319   0.017496        TCGTTT
   72322   0.029483        TTTCTG
   72323   0.039008        TTCTGT
   72324   0.031474        TCTGTG
   72326   0.026997        TGTGAG

In both of these output snippets, we see from the second column that the probability of BrdU is low (around a 1-3% chance of BrdU) so these few bases are likely from a BrdU-negative region of DNA.  In contrast, here we see the start of a read that does contain BrdU, and accordingly, the probability of BrdU at some positions is much higher:

.. code-block:: console

   >a4f36092-b4d5-47a9-813e-c22c3b477a0c chrXVI 899273 907581 fwd
   899276  0.866907        TCAAAT
   899281  0.947935        TCCACA
   899300  0.014683        TGGGAG
   899312  0.186812        TAACGG
   899320  0.934850        TTATTG

