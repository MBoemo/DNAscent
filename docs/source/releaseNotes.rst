.. _releaseNotes:

Release Notes
===============================

v2.0.0
-----------------

* Migration from HMM-based BrdU detection at every thymidine to ResNet-based detection at every thymidine,
* Significant increases to BrdU detection accuracy,
* Support for BrdU detection on GPUs,
* ``DNAscent forkSense`` to call replication origins and termination sites in both synchronously and asynchronously replicating cells at any point in S-phase,
* ``DNAscent align`` to align nanopore signals to reference,
* Significant increases to replication origin calling accuracy,
* Visualisation utility for plotting output of multiple DNAscent executables as bedgraphs.

v1.0.0
-----------------

* HMM-based BrdU detection at every thymidine,
* Improvements to BrdU detection accuracy,
* ``DNAscent train`` to train Guassian mixture models from nanopolish eventalign.

v0.1
-----------------

* HMM-based BrdU detection at ~160 thymidine-containing 6mers,
* Assignment of high- and low-BrdU regions based on Z-score, 
* Replication origin calling for early S-phase cells,
* Released with `Muller and Boemo, et al. Capturing the dynamics of genome replication on individual ultra-long nanopore sequence reads. Nature Methods 2019;16:429-436 <https://doi.org/10.1038/s41592-019-0394-y>`_.
