.. _psl:

psl
===============================

``DNAscent psl`` is a ``DNAscent`` subprogram that writes a psl file to visualise the output of ``DNAscent detect``.

Usage
-----

.. code-block:: console

   To run DNAscent psl, do:
      DNAscent psl -d /path/to/DNAscentOutput.detect -r /path/to/reference.fasta -o /path/to/psl_prefix
   Required arguments are:
     -d,--detect               path to output file from DNAscent detect,
     -r,--reference            path to genome reference in fasta format,
     -o,--output               path to output bed prefix.
   Optional arguments are:
        --threshold            probability above which a BrdU call is considered positive (default: 0.8),
        --min                  minimum read length to compute (default is 1),
        --max                  maximum read length to compute (default is Inf).

The output file from ``DNAscent detect`` should be passed using the ``-d`` flag, and the reference genome used in the alignment should be passed with the ``-r`` flag.


Output
------

The output is a psl file with each positive BrdU call marked as a tick.  These files can then be opened in IGV or the UCSC Genome Browser to visualise positive BrdU calls genome-wide.  Note that psl tracks are only plotted from the location of the first tick, so in order to visualise the portions of each read before the first BrdU call and after the last BrdU call, a placeholder tick is placed at the first and last coordinate of each read.
