.. _visualisation:

Visualisation
===============================

DNAscent supports multilevel analysis: We want users to be able to see the fork calls made by ``DNAscent forkSense`` and visualise them alongside individual base-pair resolution BrdU calls by ``DNAscent detect`` in order to see why these calls are being made.  To that end, we include a visualisation utility in ``DNAscent/utils`` that formats the output of DNAscent executables (detect and forkSense) into bedgraphs that can be visualised with IGV or the UCSC Genome Browser. You can supply this utility with the output from one or two of these executables.  If more than one is specified, the utility organises the bedgraphs so that the tracks for each read are grouped together.  

Usage
-----

.. code-block:: console

   dnascent2bedgraph.py: Converts the output of DNAscent detect and forkSense into bedgraphs.
   To run dnascent2bedgraph.py, do:
     python dnascent2bedgraph.py [arguments]
   Example:
     python dnascent2bedgraph.py -d /path/to/dnascentDetect.out -f /path/to/dnascentForksense.out -o /path/to/newBedgraphDir -n 1000 --minLength 10000
   Required arguments are at least one of the following:
     -d,--detect               path to DNAscent detect output file,
     -f,--forkSense            path to DNAscent forkSense output file.
   Required argument is:
     -o,--output               output directory which will be created.
   Optional arguments are:
        --minLength            only convert reads with specified minimum read length (in base pairs) into bedgraphs (default: 1),
        --maxLength            only convert reads with specified maximum read length (in base pairs) into bedgraphs (default: Inf),
     -n,--maxReads             maximum number of reads to convert into bedgraphs (default: Inf),
        --targets              forkSense bed file with specific reads to plot,
        --filesPerDir          maximum reads per subdirectory (default: 300).

A further example of how to use ``dnascent2bedgraph`` is given in :ref:`workflows`.

Output
------

``dnascent2bedgraph`` will create the directory you specified using the ``-o`` flag which will contain integer-numbered subdirectories.  Each of these subdirectories will contain the bedgraphs for the number of reads specified by ``--filesPerDir`` (default is 300).  If the output of more than one DNAscent executable was specified using the ``-d`` and ``-f`` flags, then the bedgraphs for each read will be grouped together so that they appear in IGV as consecutive tracks. Rather than plotting bedgraphs for every read in the sequencing run, it is sometimes useful to plot bedgraphs of just those reads that have a replication feature (e.g., an origin) called on them. When a bed file of origin calls from DNAscent forkSense is passed using the ``--targets`` flag, dnascent2bedgraph will only write bedgraphs of the reads in the bed file. 
