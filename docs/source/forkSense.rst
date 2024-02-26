.. _forkSense:

forkSense
===============================

``DNAscent forkSense`` is a ``DNAscent`` subprogram that provides a probability estimate at each thymidine that a leftward- or rightward-moving fork moved through that position during the BrdU and EdU pulses. It can assign a stall score to each fork and create replication stress signatures.

Usage
-----

.. code-block:: console

   To run DNAscent forkSense, do:
      DNAscent forkSense -d /path/to/BrdUCalls.detect -o /path/to/output.forkSense --order EdU,BrdU
   Required arguments are:
     -d,--detect               path to output file from DNAscent detect,
     -o,--output               path to output file for forkSense,
        --order                order in which the analogues were pulsed (EdU,BrdU or BrdU,EdU).
   Optional arguments are:
     -t,--threads              number of threads (default: 1 thread),
        --markAnalogues           writes analogue incorporation locations to a bed file (default: off),
        --markOrigins             writes replication origin locations to a bed file (default: off),
        --markTerminations        writes replication termination locations to a bed file (default: off),
        --markForks               writes replication fork locations to a bed file (default: off),
        --makeSignatures          writes replication stress signatures to a bed files (default: off).


The only required inputs of ``DNAscent forkSense`` is the output file produced by ``DNAscent detect`` and the order in which the analogues were pulsed.  
In the example command above, the ``--order`` flag indicates that EdU was pulsed first, and BrdU was pulsed second.  No information about the pulse length is needed.  


Output
------

Main Output File
^^^^^^^^^^^^^^^^

``DNAscent forkSense`` will produce a human-readable output file with the name and location that you specified using the ``-o`` flag.  Like the output of ``DNAscent detect``, this file starts with a short header:

.. code-block:: console

   #DetectFile /path/to/DNAscent.detect
   #Threads 1
   #Compute CPU
   #SystemStartTime 10/02/2024 13:04:33
   #Software /path/to/DNAscent
   #Version 4.0.1
   #Commit b9598a9e5bfa5f8314f92ba0f4fed39be1aee0be
   #EstimatedRegionBrdU 0.559506
   #EstimatedRegionEdU 0.202767

The fields in this header are analagous to the header from ``DNAscent detect``, but it includes two additional lines with an estimate of the thymidine-to-BrdU substitution rate in BrdU-positive regions and an estimate of the thymidine-to-EdU substitution rate in EdU-positive regions. In the example above, approximately 56% of thymidines are substituted for BrdU in BrdU-positive regions.

The rest of this file has similar formatting to that of ``DNAscent detect``.  The format for the read headers is the same.  From left to right, the tab-delimited columns indicate:

* the coordinate on the reference,
* a Boolean (0 or 1) indicating whether that position is in an EdU-positive region,
* a Boolean (0 or 1) indicating whether that position is in an BrdU-positive region.

The following example output shows an example:

.. code-block:: console

   >806d1f69-1054-4b74-8356-d935a282a22e 11 1089865 1130164 fwd
   1089873 0       0
   1089874 0       0
   1089877 0       0
   1089878 0       0
   1089879 0       0
   1089880 0       0
   1089882 0       0
   1089895 0       0
   1089899 0       0

Only reads that have at least one BrdU-positive or EdU-positive segment are written to this file. Reads with no base analogue segments called on them are omitted from this file, as 0's everywhere across these reads is implied. Note that the format of this file has changed substantially from DNAscent v2. This design decision stems from a shift in the algorithm used, as well as the desire to avoid using excess disk space with redundant information.


Bed Files
^^^^^^^^^

If the ``--markOrigins`` flag is passed, ``DNAscent forkSense`` will write the genomic region between matched leftward- and rightward-moving forks to a bed file called ``origins_DNAscent_forkSense.bed`` in the working directory.  Likewise, if the ``--markTerminations`` flag is passed, the genomic region between leftward- and rightward-moving forks moving towards each other will be recorded in a bed file called ``terminations_DNAscent_forkSense.bed``. The flag ``--markAnalogues`` will create two separate bed files: one containing the genomic location of BrdU-positive segments, and another containing the genomic location of EdU-positive segments.

If the ``--markForks flag`` is passed, two bed files will be created in the working directory. The genomic location of leftward- and rightward-moving forks will be written to separate bed files called ``leftForks_DNAscent_forkSense.bed`` and ``rightForks_DNAscent_forkSense.bed``.



All output bed files have the following space-separated columns:

* chromosome name,
* 5' boundary of the origin (or terminiation site, or fork),
* 3' boundary of the origin (or terminiation site, or fork),
* readID,
* 5' boundary of the mapped read,
* 3' boundary of the mapped read,
* strand mapped to (fwd or rev),
* fork stall score (for forks only; see below).

For origins and termination sites, the “resolution” of the calls (i.e., the third column minus the second column) will depend on your experimental setup. In synchronised early S-phase cells, the genomic distance between the 5’ and 3’ boundaries likely to be small for origins and large for termination sites, as the leftward- and rightward-moving forks should be together near the origin. In asynchronous or mid/late S-phase cells, the origin calls may appear to be a “lower’’ resolution (i.e., larger differences between the 5’ and 3’ boundaries) as the forks from a single origin will have travelled some distance before the pulses. When both forks are together at an origin, the origin bed file will record the midpoint of the analogue segment for the analogue that was pulsed first.

Fork Stalling and Pausing
^^^^^^^^^^^^^^^^^^^^^^^^^
DNAscent will assign a stall score to each called fork. These scores range from 0 (most likely unimpeded fork movement) to 1 (most likely a stall or pause). The stall score of each fork is in the last (or eigth) column of the bedfile created when the ``--markForks`` is specified. No additional input is needed; if ``--markForks`` is specified, then the fork bed files will contain stall scores. There are, however, several instances where DNAscent will decline to make a call for a fork. These include cases where DNAscent cannot see the end of the fork (e.g., if the fork runs off the read or comes together with another fork in a termination site) or if there is a nearby indel in the genomic alignment in order to avoid false positives and negatives. When this occurs, DNAscent will print a negative integer instead of a stall score clarifying why no stress call was made for this particular fork. The reason corresponding to each negative integer value is detailed in the table below.

+--------+-----------------------------------+
| Code   | Description                       |
+--------+-----------------------------------+
| -1     | Fork ends in termination site     |
+--------+-----------------------------------+
| -2     | Suspected segmentation error      |
+--------+-----------------------------------+
| -3     | Fork runs off end of read         |
+--------+-----------------------------------+
| -4     | Proximal indel > 100 bp in length |
+--------+-----------------------------------+

Stress Signatures
^^^^^^^^^^^^^^^^^
In addition to a stall score assigned to each fork, ``DNAscent forkSense`` can optionally assign an 8-dimensional stress signature to each called fork. If the ``--makeSignatures`` option is specified, two additional bed files ``leftForks_DNAscent_forkSense_stressSignatures.bed`` and ``rightForks_DNAscent_forkSense_stressSignatures.bed`` will be created in the working directory. The format is largely similar to the fork bed files above, but each line also includes an 8-dimensional stress signature for the fork call in the eight rightmost space-separated columns. From left to right, the columns are:

* chromosome name,
* 5' boundary of the origin (or terminiation site, or fork),
* 3' boundary of the origin (or terminiation site, or fork),
* readID,
* 5' boundary of the mapped read,
* 3' boundary of the mapped read,
* strand mapped to (fwd or rev),
* fork track length (in bp),
* length of the first pulsed analogue segment (in bp),
* length of the second pulsed analogue segment (in bp),
* frequency of second pulsed analogue calls in the first pulsed analogue segment,
* frequency of first pulsed analogue calls in the first pulsed analogue segment,
* frequency of first pulsed analogue calls in the second pulsed analogue segment,
* frequency of second pulsed analogue calls in the second pulsed analogue segment,
* stall score.



