.. _seeBreaks:


seeBreaks
===============================

``DNAscent seeBreaks`` is a ``DNAscent`` subprogram that that determines whether there is elevated DNA breaking at replication forks.
It was named after the three Cambridge undergraduates that developed the algorithm: Katy \ **S**\herborne, \ **E**\va Zeng, and \ **E**\mma Cohen.

``DNAscent seeBreaks`` is intended for an experimental setup whereby EdU, then BrdU (or vice versa), then excess dT are pulsed sequentially such that they are incorporated into nascent DNA to form tracks of continuous analouge incorporation in Oxford Nanopore reads. 
If these analogue tracks run all the way to the end of a sequenced molecule, this could be due to chance from random library fragmentation, or it could be due to DNA breaking at replication forks.
``DNAscent seeBreaks`` uses the nonparametric distribution of read lengths along with a mathematical model to determine whether more analogue tracks run to the end of the molecule than would be expected by chance, indicating DNA breaking at replication forks.
The algorithm will automatically tune itself to account for the read length distribution, fork speed, and analogue pulse durations that you used.

Usage
-----

.. code-block:: console

   To run DNAscent seeBreaks, do:
      DNAscent seeBreaks -l /path/to/leftForks_DNAscent_forksense.bed -r /path/to/rightForks_DNAscent_forksense.bed -d /path/to/detectOutput.bam -o /path/to/output.seeBreaks
   Required arguments are:\n"
     -l,--left                 path to leftForks file from forkSense detect with `bed` extension,
     -r,--right                path to rightFork file from forkSense detect with `bed` extension,
     -d,--detect               path to output from detect with `detect` or `bam` extension,
     -o,--output               path to output file for seeBreaks.


The only inputs required for ``DNAscent seeBreaks`` are the outputs from both ``DNAscent forkSense`` and ``DNAscent detect``. 

Output
------

``DNAscent seeBreaks`` will produce a human-readable output file with the name and location that you specified using the ``-o`` flag.  
Like other ``DNAscent`` executables, this file starts with a short header:

.. code-block:: console

   #DetectFile /path/to/DNAscent.detect
   #ForkFiles /path/to/leftForks_DNAscent_forkSense.bed /path/to/rightForks_DNAscent_forkSense.bed
   #SystemStartTime 10/08/2025 13:04:33
   #Software /path/to/DNAscent
   #Version 4.1.1
   #Commit b9598a9e5bfa5f8314f92ba0f4fed39be1aee0be

The output includes the following statistics:

.. code-block:: console

   ExpectedMean 0.35
   ExpectedStdv 0.13
   ObservedMean 0.46
   ObservedStdv 0.15
   DifferenceMean 0.05
   DifferenceStdv 0.02
   95ConfidenceInterval -0.03 0.2

- ``ExpectedMean`` and ``ExpectedStdv`` indicate the expected mean and standard deviation of the fraction of reads whose analogue tracks extend to the read end due to chance.
- ``ObservedMean`` and ``ObservedStdv`` represent the estimated mean and uncertainty of the fraction of reads whose analogue tracks extend to the read end.
- ``DifferenceMean`` and ``DifferenceStdv`` describe the difference between observed and expected values, with positive values indicating more analogue tracks extending to the read end than expected by chance.
- ``95ConfidenceInterval`` specifies the 95% confidence interval for the difference between observed and expected values. If zero lies outside this interval, it suggests significant DNA breaking at replication forks.