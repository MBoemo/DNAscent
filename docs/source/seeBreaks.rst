.. _seeBreaks:


seeBreaks
===============================

``DNAscent seeBreaks`` is a ``DNAscent`` subprogram that that determines whether there is elevated DNA breaking at replication forks.
It was named after the three Cambridge undergraduates that developed the algorithm: Katy \ **S**\herborne, \ **E**\va Zeng, and \ **E**\mma Cohen.

``DNAscent seeBreaks`` is intended for an experimental setup similar to that in our `publication <https://doi.org/10.1038/s41467-025-63168-w>`_ whereby EdU, then BrdU (or vice versa) are pulsed sequentially into replicating DNA and followed with a thymidine chase.
The analogues are incorporated into nascent DNA to form tracks of continuous analouge incorporation in Oxford Nanopore reads. 
If these analogue tracks run all the way to the end of a sequenced molecule, this could be due to chance from random library fragmentation, or it could be due to DNA breaking at replication forks.
``DNAscent seeBreaks`` uses the nonparametric distribution of read lengths along with a mathematical model to determine whether more analogue tracks run to the end of the molecule than would be expected by chance, indicating elevated DNA breaking at replication forks.
The algorithm will automatically tune itself to account for the read length distribution, fork speed, and analogue pulse durations that you used.

Flow Cell Compatibility
-----------------------

As part of ``DNAscent`` v4, ``seeBreaks`` is naturally compatible with R10.4.1 data. However, unlike other ``DNAscent`` v4 subprograms, ``DNAscent seeBreaks`` is backcompatible with legacy R9.4.1 flow cells. It is therefore perfectly acceptable to use ``DNAscent seeBreaks`` (v4.1.1) on the outputs of ``DNAscent detect`` (v3.1.2) and ``DNAscent forkSense`` (v3.1.2) generated from R9.4.1 data.

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

   #DetectFile /path/to/detectOutput.bam
   #ForkFiles /path/to/leftForks_DNAscent_forkSense.bed /path/to/rightForks_DNAscent_forkSense.bed
   #SystemStartTime 29/08/2025 16:40:36
   #Software /path/to/DNAscent/
   #Version 4.1.1
   #Commit 31214a315fcd8826fce2319a5b30348b66bfdb24

The output includes the following statistics:

.. code-block:: console

   #ExpectedReadEndFraction 0.128429
   #ExpectedReadEndFraction_StdErr 0.00931262
   #ObservedReadEndFraction 0.29308
   #ObservedReadEndFraction_StdErr 0.0125915
   #Difference 0.164357
   #Difference_StdErr 0.0157217
   #95ConfidenceInterval 0.133543 0.195172

- ``ExpectedReadEndFraction`` is the estimate of how many analogue tracks should run to the end of the read by chance. This estimate is made by using the distribution of read lengths from the output of ``DNAscent detect`` and the fork speeds from ``DNAscent forkSense``. A standard error of this estimate is determined by bootstrapping and is given by ``ExpectedReadEndFraction_StdErr``.
- ``ObservedReadEndFraction`` is the estimate of how many analogue tracks ran to the end of the read made by using the output of ``DNAscent forkSense``. A standard error of this estimate is determined by bootstrapping and is given by ``ObservedReadEndFraction_StdErr``.
- ``Difference`` is the estimate of the difference between observed and expected values, with positive values indicating more analogue tracks extending to the read end than expected by chance. The standard error of the distance estimate is given by ``Difference_StdErr``.
- ``95ConfidenceInterval`` specifies the 95% confidence interval for the difference between observed and expected values. If zero lies outside this interval, it suggests significant DNA breaking at replication forks.

Taken together, we see in the example above that 1269 forks were used, with an expected read end fraction of 0.128 (±0.0093) and an observed read end fraction of 0.293 (±0.013). The difference between observed and expected is therefore 0.164 (±0.016), with a 95% confidence interval of [0.134, 0.195]. Since zero does not lie within this confidence interval, we can conclude that there is significant DNA breaking at replication forks in this example.

Each fraction from bootstrapping on Expected and Observed is given under ``>ExpectedReadEndFractions:`` and ``>ObservedReadEndFractions:``, respectively, for those who want to investigate and/or plot these distributions. See :ref:`cookbook` for an example of how to do this in Python.