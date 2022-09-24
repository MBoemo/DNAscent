.. _index:

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
   Optional arguments are:
     -o,--output               output file name (default is index.dnascent),
     -s,--sequencing-summary   path to sequencing summary file Guppy (optional but strongly recommended).

The only required input to ``DNAscent index`` is the full path to the top-level directory containing the sequencing run's fast5 files, passed using the ``-f`` flag.  This will typically be the directory created with MinKNOW during sequencing.  An additional optional argument is the full path to the ``sequencing_summary.txt`` file, specified using the ``-s`` flag.  This file is created by Guppy during basecalling and is located in the top level directory containing the Guppy-created fastq files.  While including the sequencing summary file is optional, it is strongly recommended as it will make ``DNAscent index`` run much faster. The default behaviour of ``DNAscent index`` is to place a file called ``index.dnascent`` in the working directory.  The name of this file can be overridden using the ``-o`` flag.

Output
-------

``DNAscent index`` will put a file called ``index.dnascent`` in the current working directory (note that if you used the ``-o`` flag, then the file will have the name and location that you specified).  This file will be needed as an input to ``DNAscent detect``.
