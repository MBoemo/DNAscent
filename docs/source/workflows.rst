.. _workflows:

Workflow
===============================

The following is a full DNAscent workflow, where we'll start off after Guppy has finished running (users that need help with Guppy should refer to the `Oxford Nanopore webpages <https://nanoporetech.com/nanopore-sequencing-data-analysis>`_). The recommended Guppy basecalling configuration file for v4.0.1 is ``dna_r10.4.1_e8.2_400bps_5khz_hac.cfg``.
In particular, we assume the following:

* you have a directory of R10.4.1 Oxford Nanopore fast5 reads (which may be in subdirectories) that you want to use for detection,
* these reads have been basecalled to fastq format using Guppy (available from Oxford Nanopore),
* you have a reference/genome file (in fasta format) for your reads.

Example Workflow
----------------

Download and compile DNAscent:

.. code-block:: console

   git clone --recursive https://github.com/MBoemo/DNAscent.git
   cd DNAscent
   git checkout 4.0.1
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

This should put a file called ``index.dnascent`` in the current directory.  

You can run DNAscent detect (on 10 threads, for example) by running:

.. code-block:: console

   DNAscent detect -b alignment.sorted.bam -r /full/path/to/reference.fasta -i index.dnascent -o output.detect -t 10

Alternatively, if the system has a CUDA-compatible GPU in it, we can run ``nvidia-smi`` to get an output that looks like the following:

.. code-block:: console

   Thu Aug 20 21:06:57 2020
   +-----------------------------------------------------------------------------+
   | NVIDIA-SMI 450.51.06    Driver Version: 450.51.06    CUDA Version: 11.0     |
   |-------------------------------+----------------------+----------------------+
   | GPU  Name        Persistence-M| Bus-Id        Disp.A | Volatile Uncorr. ECC |
   | Fan  Temp  Perf  Pwr:Usage/Cap|         Memory-Usage | GPU-Util  Compute M. |
   |                               |                      |               MIG M. |
   |===============================+======================+======================|
   |   0  Tesla P100-PCIE...  On   | 00000000:05:00.0 Off |                    0 |
   | N/A   41C    P0    52W / 250W |   2571MiB / 16280MiB |     43%      Default |
   |                               |                      |                  N/A |
   +-------------------------------+----------------------+----------------------+

   +-----------------------------------------------------------------------------+
   | Processes:                                                                  |
   |  GPU   GI   CI        PID   Type   Process name                  GPU Memory | 
   |        ID   ID                                                   Usage      |
   |=============================================================================|
   |    0   N/A  N/A    178943      C   ...DNAscent_dev/bin/DNAscent     2569MiB |
   +-----------------------------------------------------------------------------+

From this, we can see that the GPU's device ID is 0 (just to the left of Tesla) so we can run:

.. code-block:: console

   DNAscent detect -b alignment.sorted.bam -r /full/path/to/reference.fasta -i index.dnascent -o output.detect -t 10 --GPU 0

Note that we're assuming the CUDA libraries for the GPU have been set up properly (see :ref:`installation`). If these libraries can't be accessed, DNAscent will splash a warning saying so and default back to using CPUs.

When ``DNAscent detect`` is finished, it will should put a file called ``output.detect`` in the current directory.  At this point, we can make bedgraphs out of the ``DNAscent detect`` output (see :ref:`visualisation`) which can also be loaded into IGV or the UCSC Genome Browser.

Lastly, we can run ``DNAscent forkSense`` on the output of ``DNAscent detect`` to measure replication fork movement.  Suppose that in our experimental protocol, we pulsed BrdU first followed by EdU.  Let's run it on four threads and specify that we want it to keep track of replication origins, forks, and termination sites:

.. code-block:: console

   DNAscent forkSense -d output.detect -o output.forkSense -t 4 --markOrigins --markTerminations --markForks --order BrdU,EdU

This will make the following files: 

* origins_DNAscent_forkSense.bed (with our origin calls),
* terminations_DNAscent_forkSense.bed (with our termination calls), 
* four bed files (leftForks_DNAscent_forkSense.bed, leftForksStressed_DNAscent_forkSense.bed, rightForks_DNAscent_forkSense.bed, rightForksStressed_DNAscent_forkSense.bed) with our fork calls,
* output.forkSense. 

We can load the bed files directly into IGV to see where origins, forks, and terminiations were called in the genome.

We can visualise (see :ref:`visualisation`) output.forkSense by turning them into bedgraphs:

.. code-block:: console

   python dnascent2bedgraph.py -d output.detect -f output.forkSense -o newBedgraphDirectory

This will create a new directory called ``newBedgraphDirectory``.  By passing both a ``forkSense`` and ``detect`` file to dnascent2bedgraph.py, the utility will convert them both into bedgraphs and organise them so that for each read, we can see the bp-resolution BrdU detection output from ``DNAscent detect`` right next to the left- and rightward-moving fork probabilities from ``DNAscent forkSense``.  These bedgraphs can then be loaded into IGV or the UCSC Genome Browser. 

Perhaps, however, we are only interested in viewing reads with origin calls on them. In this case, we can use the bed file generated above (origins_DNAscent_forkSense.bed) to specify that we only want bedgraphs of reads with origin calls on them.

.. code-block:: console

   python dnascent2bedgraph.py -d output.detect -f output.forkSense -o newBedgraphDirectory --targets origins_DNAscent_forkSense.bed
   
This strategy works equally well for any of the bed files generated by DNAscent forkSense.

Barcoding
---------

The workflow for a barcoded run is very similar to the workflow above with a few minor changes. If you're using a barcoded run that you demultiplexed with Guppy, make a fastq file for each barcode and align each of them to the reference to make as many bam files as you have barcodes. Then run ``DNAscent detect`` on the bam file for each barcode. You only have to run ``DNAscent index`` once per run, and the same ``index.dnascent`` file can be passed to ``DNAscent detect`` regardless of which barcode you're working with.

