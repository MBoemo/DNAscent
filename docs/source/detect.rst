.. _detect:

Detect
===============================

``DNAscent detect`` is a ``DNAscent`` subprogram that goes through each read and, at each thymidine position, assigns a log-likelihood of thymidine vs. BrdU.

Usage
-----

.. code-block:: console

   To run DNAscent detect, do:
     ./DNAscent detect -b /path/to/alignment.bam -r /path/to/reference.fasta -i /path/to/index.dnascent -o /path/to/output.detect
   Required arguments are:
     -b,--bam                  path to alignment BAM file,
     -r,--reference            path to genome reference in fasta format,
     -i,--index                path to DNAscent index,
     -o,--output               path to output file that will be generated.
   Optional arguments are:
     -t,--threads              number of threads (default is 1 thread),
     -q,--quality              minimum mapping quality (default is 20).
     -l,--length               minimum read length in bp (default is 100).

The main input of ``DNAscent detect`` is an alignment (bam file) between the sequence fastq from Guppy and the organism's reference genome.  This bam file should be sorted using ``samtools sort`` and indexed using ``samtools index`` so that there is a .bam.bai file in the same directory as the bam file. (Please see the example in :ref:`workflows` for details on how to do this.)  The full path to the reference genome used in the alignment should be passed using the ``-r`` flag, and the index required by the ``-i`` flag is the file created using ``DNAscent index`` (see :ref:`index`).  

Optional arguments include the number of threads, specified using the ``-t`` flag.  ``DNAscent detect`` multithreads quite well by analysing a separate read on each thread, so multithreading is recommended.  It is sometimes useful to only run ``DNAscent detect`` on reads that exceed a certain mapping quality or length threshold (as measured by the subsequence of the contig that the read maps to).  In order to do this without having to filter the bam file, DNAscent provides the ``-l`` and ``-q`` flags.  Any read in the bam file with a reference length lower than the value specificed with ``-l`` or a mapping quality lower than the value specified with ``-q`` will be ignored.

Before calling BrdU in a read, ``DNAscent detect`` must first perform a fast event alignment (see https://www.biorxiv.org/content/10.1101/130633v2 for more details).  Quality control checks are performed on these alignments, and if they're not passed, then the read fails and is ignored.  Hence, the number of reads in the output file will be slightly lower than the number of input reads.  Typical failure rates are about 5-10%, although this will vary slightly depending on the run and the genome sequenced.

Output
------

``DNAscent detect`` will produce a single human-readable output file (with the name and location that you specified using the ``-o`` flag).  Note that everything in this output file orients to the reference genome in the 5' --> 3' direction.  Each read starts with a header, which is in the following format:

.. code-block:: console

   >readID contig mappingStart mappingEnd strand

Each header always begins with a ``>`` character (therefore, an easy way to count the number of reads in the file is ``grep '>' detect.out | wc -l``).  The fields are:

* ``readID`` is a unique hexadecimal string assigned to each read by the Oxford Nanopore software,
* the read mapped between ``mappingStart`` and ``mappingEnd`` on ``contig``,
* ``strand`` either takes the value ``fwd``, indicating that the read mapped to the forward strand, or ``rev`` indicating that the read mapped to the reverse complement strand.

The following shows an example header for a read that to the reverse strand between 239248 and 286543 on chrII.

.. code-block:: console

   >c602f23f-e892-42ba-8140-da949abafbdd chrII 239248 286543 rev

Below the header, each line corresponds to the position of a thymidine in the read.  There are four tab-separated columns:

* the coordinate on the reference,
* log-likelihood that the thymidine is actually BrdU (positive values look more like BrdU, negative values look more like thymidine),
* 6mer on the reference,
* Guppy-basecalled 6mer on the query (or its reverse complement if the read mapped to the reverse strand).


Consider the following example:

.. code-block:: console

   >c602f23f-e892-42ba-8140-da949abafbdd chrII 239248 286543 rev
   239265  -3.695165       TATGCA  TATGCA
   239267  -2.448894       TGCAGA  TGCAGA
   239269  -2.806177       CAGATA  CAGATA
   239270  -2.734358       AGATAA  AGATAA

Here, we're looking at the sequence TATGCAGATAA on the reference genome.  Because this read maps to the reverse complement, a call is made at every A (instead of T) on the reference.  If instead we looked at a read that mapped to the forward strand, an example would be:

.. code-block:: console

   >f19253eb-a773-439d-a824-47ceea355109 chrIV 12407 56248 fwd
   12433   -10.100426      TTCTTT  TTCTTT
   12434   -10.052507      TCTTTC  TCTTTC
   12436   -10.244013      TTTCCA  TTTCCA
   12437   -11.881122      TTCCAT  TTCCAT
   12438   -12.543487      TCCATG  TCCATG
