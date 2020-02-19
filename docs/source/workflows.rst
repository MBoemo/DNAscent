.. _workflows:

Sample Workflow
===============================

The following is a full DNAscent workflow, where we'll start off after Guppy has finished running (users that need help with Guppy should refer to the `Oxford Nanopore webpages <https://nanoporetech.com/nanopore-sequencing-data-analysis>`_).  In particular, we assume the following:

* you have a directory of 1D R9.5 450bp/s Oxford Nanopore fast5 reads (which may be in subdirectories) that you want to use for detection,
* these reads have been basecalled to fastq format using Albacore or Guppy (available from Oxford Nanopore),
* you have a reference/genome file (in fasta format) for your reads.

Download and compile DNAscent:

.. code-block:: console

   git clone --recursive https://github.com/MBoemo/DNAscent.git
   cd DNAscent
   make
   cd ..

Concatenate the fastq files from Guppy:

.. code-block:: console

   cat /path/to/GuppyOutDirectory/*.fastq > reads.fastq

Align the reads with `minimap2 <https://github.com/lh3/minimap2>`_ and sort with `samtools <http://www.htslib.org/>`_:

.. code-block:: console

   minimap2 -ax map-ont -o alignment.sam /path/to/reference.fasta reads.fastq
   samtools view -Sb -o alignment.bam alignment.sam
   samtools sort alignment.bam alignment.sorted
   samtools index alignment.sorted.bam

Now we're ready to use DNAscent.  Let's index the run:

.. code-block:: console

   DNAscent index -f /full/path/to/fast5 -s /full/path/to/GuppyOutDirectory/sequencing_summary.txt

This should put a file called ``index.dnascent`` in the current directory.  You can run DNAscent detect (on 10 threads, for example) by running:

.. code-block:: console

   DNAscent detect -b alignment.sorted.bam -r /full/path/to/reference.fasta -i index.dnascent -o output.detect -t 10

This should put a file called ``output.detect`` in the current directory.  We can look at the individual positive BrdU calls by running:

.. code-block:: console

   DNAscent psl -d output.detect -r /full/path/to/reference.fasta -o output

The resulting file ``output.psl`` can be loaded into IGV or the UCSC Genome Browser.

We can measure the regions of each read with high BrdU incorporation by using ``DNAscent regions``:

.. code-block:: console

   DNAscent regions -d output.detect -o output.regions

Lastly, we can take a look at the DNAscent regions results in IVG by generating a bedgraph file for each read.  Note that the bedgraph format requires that we make a separate file for each read, so we'll want to make a separate directory for this.

.. code-block:: console

   mkdir bedgraphs
   cd bedgraphs
   python /path/to/DNAscent/scripts/regions2bedgraph.py ../output.regions

To make things a little more manageable, regions2bedgraph.py puts the bedgraph files into numbered subdirectories such that each contains about 400 reads. These bedgraph files can then be viewed in IGV or the UCSC Genome Brower.
