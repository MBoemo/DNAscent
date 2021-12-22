.. _index:

Index
===============================

``DNAscent index`` is a ``DNAscent`` subprogram that creates a map between Oxford Nanopore readIDs and fast5 files.  This allows ``DNAscent detect`` to scan through bam files and pull out the relevant signal information for each read.

Usage
-----

.. code-block:: console

   To run DNAscent index, do:
      DNAscent index -f /path/to/fast5Directory
   Required arguments are:
     -f,--files                path to fast5 files,
     -s,--sequencing-summary   path to sequencing summary file Guppy.
   Optional arguments are:
     -o,--output               output file name (default is index.dnascent),
        --GridION              account for the different sequencing summary format used by in-built GridION basecalling.

The first required input to ``DNAscent index`` is the full path to the top-level directory containing the sequencing run's fast5 files, passed using the ``-f`` flag.  This will typically be the directory created with MinKNOW during sequencing.  The second required input is the full path to the ``sequencing_summary.txt`` file, specified using the ``-s`` flag.  This file is created by Guppy during basecalling and is located in the top level directory containing the Guppy-created fastq files.  (Note that as of v2.0.3, providing the sequencing summary file is now required where it was optional in previous releases.) The default behaviour of ``DNAscent index`` is to place a file called ``index.dnascent`` in the working directory.  The name of this file can be overridden using the ``-o`` flag. Note that the in-built version of Guppy on the GridION produces a sequencing summary file with a slightly different format than the version of Guppy available on the ONT Community website. Users can correct for this change of format from the GridION by adding the ``--GridION`` flag.

Output
-------

``DNAscent index`` will put a file called ``index.dnascent`` in the current working directory (note that if you used the ``-o`` flag, then the file will have the name and location that you specified).  This file will be needed as an input to ``DNAscent detect``.
