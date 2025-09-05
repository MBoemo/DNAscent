.. DNAscent documentation master file, created by
   sphinx-quickstart on Mon Feb 26 13:46:59 2024.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

DNAscent
====================================

.. toctree::
   :maxdepth: 1
   :caption: Contents:

   installation
   index_exe
   detect
   forkSense
   seeBreaks
   visualisation
   workflows
   cookbook
   releaseNotes


Overview
--------

DNAscent is software designed to detect the base analogues BrdU and EdU in Oxford Nanopore reads.  In an experimental setup where BrdU and EdU are incorporated into nascent DNA by replication forks, this software can be used to answer questions that were traditionally answered by DNA fibre analysis.  DNAscent can also call the genomic positions of stalled and stressed replication forks for use as a replication stress assay.

DNAscent v4.1.1 supports sequencing data collected on Oxford Nanopore R10.4.1 flow cells. Users wishing to analyse data acquired on legacy R9.4.1 flow cells should roll back to DNAscent v3.1.2 as v4.1.1 is not back-compatible with R9.4.1 flow cells. As R9.4.1 flow cells have been deprecated by Oxford Nanopore,
previous versions of DNAscent designed for R9.4.1 flow cells (v3.1.2 and below) are no longer under active development.

The Oxford Nanopore Flongle, MinION, GridION, and PromethION platforms are all supported.

DNAscent is under active development by the `Boemo Group <https://www.boemogroup.org/>`_ based in the `Department of Pathology, University of Cambridge <https://www.path.cam.ac.uk/>`_.  We aim to push regular updates and improvements and incorporating new functionality is an active area of our computational research.


Developer Recommendations
-------------------------

* **R9.4.1 or R10.4.1?** R10.4.1. Users wanting to use legacy R9.4.1 data should roll back to DNAscent v3.1.2. All future development will be for the R10.4.1 pore.
* **FAST5 or POD5?** POD5. Dorado and DNAscent both perform better with POD5 compared to FAST5.
* **Guppy or Dorado?** Basecalling with Dorado to bam format is our recommended workflow. It is faster than Guppy, has the newest basecalling models, and writes informative information to the output bam file which DNAscent can make use of.
* **Singularity image or compile from source?** Singularity image. We load it with DNAscent's dependencies so that all you need is a valid NVIDIA driver. 
* **Human-readable or modbam output format from DNAscent detect?** We recommend modbam for most users. It produces about a 6x savings in drive space, can be visualised directly by genomics viewers, and makes some downstream analyses easier. 
* **GPU or CPU?** GPU for DNAscent detect (along with as many CPUs as your GPU node has available) and CPU for DNAscent align and forkSense. DNAscent index is quick and should run in a few seconds on a single CPU.

Publications
------------

If you use DNAscent for your research, please cite the publication appropriate for the version you used:

**v3 & v4:** Jones MJK,  Rai SK,  Pfuderer PL, Bonfim-Melo A, Pagan JK, Clarke PR, Totanes FIG, Merrick CJ, McClelland SE, Boemo MA. A high-resolution, nanopore-based artificial intelligence assay for DNA replication stress in human cancer cells. [`Journal DOI <https://doi.org/10.1038/s41467-025-63168-w>`_]

**v3:** Totanes FIG,  Gockel J,  Chapman SE, Bartfai R, Boemo MA, Merrick CJ. Replication origin mapping in the malaria parasite Plasmodium falciparum. [`Journal DOI <https://doi.org/10.1093/nar/gkad093>`_]

**v2:** Boemo, MA DNAscent v2: Detecting replication forks in nanopore sequencing data with deep learning. BMC Genomics 2021;22:430. [`Journal DOI <https://doi.org/10.1186/s12864-021-07736-6>`_]

**v1:** Muller CA, Boemo MA, Spingardi P, Kessler, BM, Kriaucionis S, Simpson JT, Nieduszynski CA. Capturing the dynamics of genome replication on individual ultra-long nanopore sequence reads. Nature Methods 2019;16:429-436. [`Journal DOI <https://doi.org/10.1038/s41592-019-0394-y>`_]

Bugs, Questions, and Comments
-----------------------------

Should any bugs arise or if you have any questions about usage, please raise a `GitHub issue <https://github.com/MBoemo/DNAscent/issues>`_. Your feedback is an important part of the development process. If you have comments or suggestions to improve the software or the documentation, please Email Michael Boemo at mb915@cam.ac.uk.
