.. _index_exe:

index
===============================

``DNAscent index`` is a ``DNAscent`` subprogram that creates a map between Oxford Nanopore readIDs and FAST5/POD5 files.  This allows ``DNAscent detect`` to scan through bam files and pull out the relevant signal information for each read.

Usage
-----

.. code-block:: console

   To run DNAscent index, do:
      DNAscent index -f /path/to/directory
   Required arguments are:\n"
      -f,--files                full path to fast5 or pod5 files.
   Optional arguments are:\n"
      -s,--sequencing-summary   (legacy) path to sequencing summary file from using Guppy on fast5 files,
      -o,--output               output file name (default is index.dnascent).

The one required input to ``DNAscent index`` is the full path to the top-level directory containing the sequencing run's FAST5 files or POD5 files (passed using the ``-f`` flag). It is permissible to pass a directory containing both FAST5 and POD5 files. We recommend passing data to DNAscent in POD5 format. 

For users that are using legacy FAST5 format and Guppy, you can either:

* Pass only the FAST5 directory using the ``-f`` flag. DNAscent will iterate through each FAST5 file to build the index (expected to be slow).
* Pass the FAST5 directory using the ``-f`` flag and the ``sequencing_summary.txt`` file from Guppy using the ``-s`` flag. This will be much faster.  

The default behaviour of ``DNAscent index`` is to place a file called ``index.dnascent`` in the working directory. The name of this file can be overridden using the ``-o`` flag.

Output
-------

``DNAscent index`` will put a file called ``index.dnascent`` in the current working directory (note that if you used the ``-o`` flag, then the file will have the name and location that you specified).  This file will be needed as an input to ``DNAscent detect``.
