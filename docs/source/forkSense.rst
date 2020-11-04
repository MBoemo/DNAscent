.. _forkSense:

forkSense
===============================

``DNAscent forkSense`` is a ``DNAscent`` subprogram that provides a probability estimate at each thymidine that a leftward- or rightward-moving fork moved through that position during the BrdU pulse.

Usage
-----

.. code-block:: console

   To run DNAscent forkSense, do:
      DNAscent forkSense -d /path/to/BrdUCalls.detect -o /path/to/output.forkSense
   Required arguments are:
     -d,--detect               path to output file from DNAscent detect,
     -o,--output               path to output file for forkSense.
   Optional arguments are:
     -t,--threads              number of threads (default: 1 thread),
        --markOrigins             writes replication origin locations to a bed file (default: off),
        --markTerminations        writes replication termination locations to a bed file (default: off),
        --markForks               writes replication fork locations to a bed file (default: off).


The only required input of ``DNAscent forkSense`` is the output file produced by ``DNAscent detect``.  Note that the detect file must have been produced using the v2.0 ResNet algorithm; ``DNAscent forkSense`` is not compatible with legacy HMM-based detection.

If the ``--markOrigins`` flag is passed, ``DNAscent forkSense`` will use detected leftward- and rightward-moving forks to infer the locations of fired replication origins and write these to a bed file called ``origins_DNAscent_forkSense.bed`` in the working directory.  Likewise, if the ``--markTerminations`` flag is passed, termination sites will be recorded in a bed file called ``terminations_DNAscent_forkSense.bed``.

Output
------

If ``--markOrigins`` and/or ``--markTerminations`` were used, the resulting bed files has one called origin (for origins_DNAscent_forkSense.bed) or termination site (for terminations_DNAscent_forkSense.bed) per line and, in accordance with bed format, have the following space-separated columns:

* chromosome name,
* 5' boundary of the origin (or terminiation site),
* 3' boundary of the origin (or terminiation site),
* read header of the read that the call came from (similar to those in the output file of ``DNAscent detect``).

Note that the "resolution" of the calls (i.e., the third column minus the second column) will depend on your experimental setup.  In synchronised early S-phase cells, this difference for origin calls is likely to be small as the leftward- and rightward-moving forks from a fired origin are nearby one another.  In asynchronous or mid/late S-phase cells, the difference is likely to be larger as the forks from a single origin will have travelled some distance before the BrdU pulse.  The bed files only specify the region between matching leftward- and rightward-moving forks.  Any subsequent assumptions (such as assuming uniform fork speed and placing the origin in the middle of that region) are left to the user.

The output of ``DNAscent forkSense`` is a file with similar formatting to that of ``DNAscent detect``.  The format for the read headers is the same.  From left to right, the tab-delimited columns indicate:

* the coordinate on the reference,
* probability that a leftward-moving fork passed through that coordinate during a BrdU pulse,
* probability that a rightward-moving fork passed through that coordinate during a BrdU pulse.

A low probability in both the second and third columns suggests it was unlikely that a fork passed through that position during the pulse.

The following example output shows the end of a read that was passed through by a leftward-moving fork:

.. code-block:: console

   >22c8a674-ed0e-475f-9c54-cb185299d923 chrII 173332 210452 fwd
   173339  0.687217        0.062620
   173341  0.687217        0.062620
   173342  0.687217        0.062620
   173345  0.687217        0.062620
   173347  0.687217        0.062620
   173348  0.687217        0.062620
   173349  0.743986        0.045767
   173358  0.743986        0.045767
   173375  0.743986        0.045767
   173377  0.743986        0.045767
   173378  0.743986        0.045767
   173381  0.743986        0.045767
   173382  0.806924        0.038138
   173383  0.806924        0.038138
   173387  0.806924        0.038138
   173390  0.806924        0.038138
   173392  0.806924        0.038138
   173393  0.806924        0.038138
   173398  0.846875        0.032027
   173402  0.846875        0.032027
   173404  0.846875        0.032027
   173406  0.846875        0.032027
   173407  0.846875        0.032027
   173417  0.846875        0.032027
   173418  0.906748        0.028587
   173419  0.906748        0.028587
   173423  0.906748        0.028587
   173425  0.906748        0.028587
   173426  0.906748        0.028587
   173428  0.906748        0.028587
   173441  0.909755        0.029341
   173445  0.909755        0.029341
   173446  0.909755        0.029341
   173449  0.909755        0.029341
   173450  0.909755        0.029341
   173451  0.909755        0.029341
   173454  0.907803        0.029983


