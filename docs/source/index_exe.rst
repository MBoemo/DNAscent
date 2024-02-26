.. _index_exe:

index
===============================

``DNAscent index`` is a ``DNAscent`` subprogram that creates a map between Oxford Nanopore readIDs and fast5 files.  This allows ``DNAscent detect`` to scan through bam files and pull out the relevant signal information for each read.

Usage
-----

.. code-block:: console

   To run DNAscent index, do:
      DNAscent index -f /path/to/fast5Directory
   Required arguments are:
     -f,--files                path to fast5 files.
     -s,--sequencing-summary   path to sequencing summary file Guppy.
   Optional arguments are:
     -o,--output               output file name (default is index.dnascent),
        --GridION              account for the different sequencing summary format used by in-built GridION basecalling.

The required inputs to ``DNAscent index`` are the full path to the top-level directory containing the sequencing run's fast5 files (passed using the ``-f`` flag) and the path to the ``sequencing_summary.txt`` file (specified using the ``-s`` flag).  
``sequencing_summary.txt`` is created by Guppy during basecalling and is located in the top level directory containing the Guppy-created fastq files.  
The default behaviour of ``DNAscent index`` is to place a file called ``index.dnascent`` in the working directory.  The name of this file can be overridden using the ``-o`` flag.

Output
-------

``DNAscent index`` will put a file called ``index.dnascent`` in the current working directory (note that if you used the ``-o`` flag, then the file will have the name and location that you specified).  This file will be needed as an input to ``DNAscent detect``.
