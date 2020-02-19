.. DNAscent documentation master file, created by
   sphinx-quickstart on Fri Feb  7 18:58:49 2020.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

DNAscent
====================================

.. toctree::
   :maxdepth: 1
   :caption: Contents:

   installation
   index
   detect
   regions
   psl
   workflows

Overview
--------

DNAscent is software designed to detect the modified base BrdU in Oxford Nanopore reads.  It uses a hidden Markov approach to differentiate between BrdU and thymidine using the raw nanopore signal.  In an experimental setup where BrdU is incorporated into nascent DNA by replication forks, this software can be used to answer questions that were traditionally answered by DNA fibre analysis.  

DNAscent is under active development by the `Boemo Group <https://www.boemogroup.org/>`_ based in the `Department of Pathology, University of Cambridge <https://www.path.cam.ac.uk/>`_.  We aim to push regular updates and improvements, and incorporating new functionality is an active area of our computational research.


Publications
------------

Please cite the following publication if you use DNAscent for your research: 

Muller CA, Boemo MA, Spingardi P, Kessler, BM, Kriaucionis S, Simpson JT, Nieduszynski CA. Capturing the dynamics of genome replication on individual ultra-long nanopore sequence reads. Nature Methods 2019;16:429-436. [`Journal DOI <https://doi.org/10.1038/s41592-019-0394-y>`_]

Bugs, Questions, and Comments
-----------------------------

Should any bugs arise or if you have any questions about usage, please raise a `GitHub issue <https://github.com/MBoemo/DNAscent/issues>`_. If you have comments or suggestions to improve the software or the documentation, please Email Michael Boemo at mb915@cam.ac.uk.
