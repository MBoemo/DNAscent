.. _releaseNotes:

Release Notes
===============================

v4.0.1
-----------------

This is a pre-release. It is possible (maybe even probable) that the deep learning model and other parts of the software will be re-tuned in the coming weeks as we listen to community feedback. 

* All DNAscent executables now work on Oxford Nanopore R10.4.1 pores. ``DNAscent detect`` detects the thymidine analogues BrdU and EdU in single molecules as in v3.0.2 and v3.1.2.
* POD5 is not yet supported but POD5 support is planned for v4.0.2 which will be the LTS release. POD5 can be converted to FAST5 using the ``pod5 convert to_fast5`` utility (https://pypi.org/project/pod5/#pod5-convert-to_fast5).
* DNAscent's deep learning models were trained on R10.4.1 flow cells with 5kHz sampling which is now the default sampling rate on Oxford Nanopore platforms. Using DNAscent with reads sequenced with the older 4kHz sampling rate is not recommended.
* ``DNAscent index`` still uses the sequencing_summary.txt file from Guppy and legacy versions of Guppy are available on the Oxford Nanopore Community webpage. Compatibility with Dorado is planned for v4.0.2 LTS.
* Tensorflow updated to 2.12.0 and, correspondingly, GPU usage now requires CUDA 11.8 and cuDNN 8.9.
* Useage is otherwise identical to v3.1.2.
* Training data for this release was provided by Mathew Jones at the University of Queensland. This software was developed in collaboration with the `Jones Lab <https://researchers.uq.edu.au/researcher/25051>`_, `Merrick Lab <https://www.path.cam.ac.uk/directory/catherine-merrick>`_, and `McClelland Lab <https://www.bartscancer.london/staff/professor-sarah-mcclelland/>`_ to whom we are grateful for their collaboration and support. We are particularly grateful to them for supporting the release of the software to the community ahead of publication.

v3.1.2
-----------------

* ``DNAscent forkSense`` now assigns a stall score to each called fork,
* ``DNAscent forkSense`` can assign a replication stress signature to each called fork,
* ``DNAscent detect`` no longer outputs the reference 6mer corresponding to each thymidine position in order to reduce output file size,
* improvements to fork calling and segmentation,
* Released with `Jones MJK,  Rai SK,  Pfuderer PL, Bonfim-Melo A, Pagan JK, Clarke PR, McClelland SE, Boemo MA. A high-resolution, nanopore-based artificial intelligence assay for DNA replication stress in human cancer cells. bioRxiv <https://doi.org/10.1101/2022.09.22.509021>`_.

v3.0.2
-----------------

* ``DNAscent detect`` now detects two different thymidine analogues, BrdU and EdU, in the same molecule,
* ``DNAscent forkSense`` now uses the spatial patterning of EdU and BrdU to determine fork direction as in DNA fibre,
* dnascent2bedgraph utility updated to plot both EdU and BrdU tracks in genome browsers,
* ``DNAscent regions`` is now deprecated and has been fully superceded by ``DNAscent forkSense``,
* ``DNAscent psl`` is now deprecated as reads can be more comprehensively plotted using the dnascent2bedgraph utility,
* Migration from Tensorflow 1.14 to 2.4.1 and, correspondingly, GPU usage now requires CUDA 11 and cuDNN 8,
* Released with `Totanes FIG,  Gockel J,  Chapman SE, Bartfai R, Boemo MA, Merrick CJ. Replication origin mapping in the malaria parasite Plasmodium falciparum. bioRxiv <https://doi.org/10.1101/2022.07.27.501677>`_.

v2.0.0
-----------------

* Migration from HMM-based BrdU detection at every thymidine to ResNet-based detection at every thymidine,
* Significant increases to BrdU detection accuracy,
* Support for BrdU detection on GPUs,
* ``DNAscent forkSense`` to call replication origins and termination sites in both synchronously and asynchronously replicating cells at any point in S-phase,
* ``DNAscent align`` to align nanopore signals to reference,
* Significant increases to replication origin calling accuracy,
* Visualisation utility for plotting output of multiple DNAscent executables as bedgraphs,
* Released with `Boemo, MA. DNAscent v2: Detecting replication forks in nanopore sequencing data with deep learning. BMC Genomics 2021;22:430 <https://doi.org/10.1186/s12864-021-07736-6>`_.

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
