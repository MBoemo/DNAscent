.. _refSignal:

refSignal
===============================

``DNAscent refSignal`` is a ``DNAscent`` subprogram that builds a genome-wide reference nanopore signal profile from a sorted alignment BAM file.  For every reference position that has at least one read aligned to it, it estimates the expected pore current as a Gaussian distribution parameterised by a mean and standard deviation.  This profile can subsequently be used by ``DNAscent refScan`` to identify positions where a read's signal deviates from the baseline, which is consistent with the presence of a base modification.

Usage
-----

.. code-block:: console

   To run DNAscent refSignal, do:
      DNAscent refSignal -b /path/to/sorted.bam -r /path/to/reference.fasta -i /path/to/index.dnascent -o /path/to/output
   Required arguments are:
     -b,--bam                  path to sorted alignment BAM file,
     -r,--reference            path to genome reference in fasta format,
     -i,--index                path to DNAscent index,
     -o,--output               path to output file.
   Optional arguments are:
     -t,--threads              number of threads (default is 1 thread),
     -q,--quality              minimum mapping quality (default is 20),
     -l,--length               minimum read length in bp (default is 1000),
     -d,--maxDepth             maximum number of reads used per position to estimate
                               the Gaussian (default is 30; set to 0 for no limit).

The main input of ``DNAscent refSignal`` is a sorted BAM file produced by aligning basecalled reads to a reference genome.  The reference genome passed with ``-r`` should be the same FASTA file used to produce that alignment.  The index required by ``-i`` is the file created using ``DNAscent index`` (see :ref:`index_exe`).

``DNAscent refSignal`` is intended to be run on a dataset of reads that are free of, or have a negligible level of, base modifications.  The resulting signal profile then serves as the unmodified baseline against which experimental reads can be compared using ``DNAscent refScan``.  A typical workflow would therefore be:

1. Sequence an unmodified control sample and align the reads to the reference.
2. Run ``DNAscent refSignal`` on the control BAM to produce the reference signal profile.
3. Run ``DNAscent refScan`` on the experimental BAM, providing the profile from step 2.

The number of threads is specified using the ``-t`` flag.  ``DNAscent refSignal`` parallelises by processing a batch of reads simultaneously, so multithreading is recommended.

For high-coverage datasets, the ``-d`` flag caps the number of reads that contribute to the Gaussian estimate at any single reference position.  Once a position has accumulated ``-d`` observations it is considered saturated and further reads covering it are skipped early, before signal fetching and event alignment, so they incur no meaningful computational cost.  The default of 30 is sufficient for a reliable mean and standard deviation; set ``-d 0`` to apply no limit.

Method
------

For each read that passes the quality filters, ``DNAscent refSignal`` fetches the raw pore current signal from the corresponding FAST5 or POD5 file and performs a signal normalisation and event alignment identical to the one used by ``DNAscent align``.  The event alignment maps each segment of the raw signal to the reference position it most likely corresponds to.

For every reference position that receives at least one aligned signal observation, ``DNAscent refSignal`` accumulates statistics using Welford's online algorithm.  Each read contributes exactly one data point per covered position,the mean of the raw signal values aligned to that position, so that every read is weighted equally regardless of dwell time.  After all reads have been processed, the sample mean and sample standard deviation of these per-read means are written to the output file.

Output
------

``DNAscent refSignal`` writes a single tab-separated output file to the path specified with ``-o``.  The file begins with a short header in which every line starts with a hash (``#``) character:

.. code-block:: console

   #Alignment /path/to/sorted.bam
   #Genome /path/to/reference.fasta
   #Index /path/to/index.dnascent
   #Threads 4
   #MappingQuality 20
   #MappingLength 1000
   #MaxDepth 30
   #Version 4.2.2
   #Commit 4cf80a7b89bdf510a91b54572f8f94d3daf9b167

You can access the header with ``grep '#' /path/to/output``.

Below the header is one data row per covered reference position, sorted by contig name and then by position.  The columns are:

* ``contig`` — the name of the reference contig,
* ``position`` — the 0-based coordinate on that contig,
* ``mean_signal`` — the mean pore current (in normalised pA) across all reads covering this position,
* ``std_signal`` — the sample standard deviation of the per-read signal means at this position,
* ``coverage`` — the number of reads that contributed a signal observation at this position.

An example excerpt is:

.. code-block:: console

   chr1    10042     101.342       5.217       38
   chr1    10043     98.714        4.883       38
   chr1    10044     103.561       6.104       37
   chr1    10045     99.027        5.640       39

Positions with ``std_signal`` equal to zero were covered by exactly one read, and their standard deviation cannot be estimated.  ``DNAscent refScan`` skips such positions when computing window scores.

A log file is also written alongside the output.  Its name is derived from the output filename with a ``.refSignal.log`` suffix.  Any reads that could not be located in the DNAscent index are reported there.
