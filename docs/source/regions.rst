.. _regions:

Regions
===============================

``DNAscent regions`` is a ``DNAscent`` subprogram that interprets the output of ``DNAscent detect`` to call regions of high and low BrdU incorporation.

Note that as of v2.0, ``DNAscent regions`` has been largely superceded by ``DNAscent forkSense`` and the increased accuracy of ``DNAscent detect`` makes visualising BrdU incorporation in regions mostly unnecessary.  However, it is still included to avoid breaking legacy workflows, and it does still have some uses as explained below.

Usage
-----

.. code-block:: console

   To run DNAscent regions, do:
      DNAscent regions -d /path/to/DNAscentOutput.detect -o /path/to/DNAscentOutput.regions
   Required arguments are:
     -d,--detect               path to output file from DNAscent detect,
     -o,--output               path to output directory for bedgraph files.
   Optional arguments (if used with default ResNet-based detect) are:
     -r,--resolution           number of thymidines in a region (default is 10).
   Optional arguments (if used with HMM-based detect) are:
        --threshold            probability above which a BrdU call is considered positive (default: 0.8),
     -c,--cooldown             minimum gap between positive analogue calls (default: 4),
     -r,--resolution           minimum length of regions (default is 100 bp),
     -p,--probability          override probability that a thymidine 6mer contains a BrdU (default: automatically calculated),
     -z,--zScore               override zScore threshold for BrdU call (default: automatically calculated).

The only required input of ``DNAscent regions`` is the output file produced by ``DNAscent detect``. ``DNAscent regions`` will first look through the detect file and determine the approximate fraction of thymidines replaced by BrdU in BrdU-positive regions.  Using this probability, a z-score is assigned to each window (100 bp wide by default, but this can be changed using the ``-r`` flag) to indicate whether there is more or less BrdU than would be expected for an average BrdU-positive region.  Naturally, some regions will be BrdU-positive but will have a substitution rate lower than average for BrdU-positive regions. Hence, ``DNAscent regions`` determines an appropriate boundary threshold between BrdU-positive regions and thymidine regions and rescales all of the z-scores so that this boundary is 0. ``DNAscent regions`` will calculate these values for you, but they can be overridden with the  ``-p`` and ``-z`` flags, though this is generally not recommended.  The exceptions are runs with 0% BrdU or runs where a high BrdU incorporation is expected along the entirety of each read. This is because these parameters are computed assuming that there are two populations (BrdU-positive and thymidine-only segments of DNA).

In order to determine regions of high and low BrdU incorporation, ``DNAscent regions`` needs to count positive BrdU calls.  By default, a thymidine is considered to be BrdU if it was scored with a probability higher than 0.8 by ``DNAscent detect``.  This value was tuned in-house to optimise signal-to-noise, but it can be changed with the ``--threshold`` flag.  Likewise, some care has to be given to how positive calls are counted, as BrdU can sometimes shift the signal of neighbouring thymidines.  To prevent artefacts from overcounting while minimising undercounting, the default behaviour is to only make a positive call at most every 4 bases, though this can be changed with the ``-c`` flag.


Output
------

The output of DNAscent regions is a file with similar formatting to that of ``DNAscent detect``.  The format for the read headers is the same.  From left to right, the tab-delimited columns indicate:

* the start of the region,
* the end of the region,
* the z-score,
* the string "BrdU" if the score is positive and "Thym" if the score is negative.

A large positive z-score indicates high BrdU incorporation in that region, and a large negative score indicates very little BrdU incorporation in that region.  An example output is as follows:

.. code-block:: console

   >bfdc06e0-001f-41f7-bbea-f2f6785a3860 chrI 0 28066 fwd
   62      167     -2.38086        Thym
   173     276     -2.38086        Thym
   283     388     -2.27466        Thym
   393     499     -2.00741        Thym
   501     605     -2.00741        Thym
   606     708     -2.48397        Thym
   713     817     -2.27466        Thym

Note that the region width may sometimes vary slightly from the value specified. The region width is designated as the coordinate of the first thymidine greater than the window width (100 bp by default) from the starting coordinate.  In order to guard against assigning a score to regions with very few thymidines, ``DNAscent regions`` will also extend the region until at least 10 calls are considered.
