# DNAscent
[![Documentation Status](https://readthedocs.org/projects/dnascent/badge/?version=latest)](https://dnascent.readthedocs.io/en/latest/?badge=latest)

DNAscent is software designed to detect the base analogues BrdU and EdU in single molecules of DNA sequenced on the Oxford Nanopore platform. In an experimental setup where BrdU and EdU are incorporated into nascent DNA by replication forks, this software can be used to answer questions that were traditionally answered by DNA fibre analysis. DNAscent can also call the genomic positions of stalled and stressed replication forks for use as a replication stress assay.

DNAscent v4.0.1 supports sequencing data collected on Oxford Nanopore R10.4.1 flow cells. The Oxford Nanopore Flongle, MinION, GridION, and PromethION platforms are all supported.

DNAscent is under active development by the [Boemo Group](https://www.boemogroup.org/) based in the [Department of Pathology, University of Cambridge](https://www.path.cam.ac.uk/).

## Downloading and Compiling DNAscent
Clone the DNAscent repository with the recursive flag so that the dependencies are cloned as well.
```shell
git clone --recursive https://github.com/MBoemo/DNAscent.git
```
The DNAscent directory will appear in your current directory.  Switch to the latest tagged version and compile the software by running:
```shell
cd DNAscent
git checkout 4.0.1
make
```
This will put the DNAscent executable into the DNAscent/bin directory.  A typical compile time for DNAscent and its dependencies is 5 minutes.

## Documentation
Please see the [documentation](https://dnascent.readthedocs.io) for detailed usage instructions, descriptions of DNAscent's subprograms, and an example workflow.

## Citation
Please cite our publications if you use DNAscent for your research:
- Jones MJK,  Rai SK,  Pfuderer PL, Bonfim-Melo A, Pagan JK, Clarke PR, McClelland SE, Boemo MA. A high-resolution, nanopore-based artificial intelligence assay for DNA replication stress in human cancer cells. bioRxiv. [[bioRxiv](https://doi.org/10.1101/2022.09.22.509021)]
- Totanes FIG,  Gockel J,  Chapman SE, Bartfai R, Boemo MA, Merrick CJ. Replication origin mapping in the malaria parasite Plasmodium falciparum. bioRxiv. [[bioRxiv](https://doi.org/10.1101/2022.07.27.501677)]
- Boemo, MA. DNAscent v2: Detecting replication forks in nanopore sequencing data with deep learning. *BMC Genomics* 2021;22:430. [[Journal Link](https://doi.org/10.1186/s12864-021-07736-6)]
- Muller CA, Boemo MA, Spingardi P, Kessler BM, Kriaucionis S, Simpson JT, Nieduszynski CA. Capturing the dynamics of genome replication on individual ultra-long nanopore sequence reads. *Nature Methods* 2019;16:429-436. [[Journal Link](https://www.nature.com/articles/s41592-019-0394-y)]

## Questions and Bugs
Should any bugs arise or if you have basic usage questions, please raise a [GitHub issue](https://github.com/MBoemo/DNAscent/issues). For more detailed discussions or collaborations, please Email Michael Boemo at mb915@cam.ac.uk.

## Legacy Flow Cells
Users wishing to analyse data acquired on legacy R9.4.1 flow cells should roll back to DNAscent v3.1.2 as v4.0.1 and above are not back-compatible with R9.4.1 flow cells. As R9.4.1 flow cells have been deprecated by Oxford Nanopore,
previous versions of DNAscent designed for R9.4.1 flow cells (v3.1.2 and below) are no longer under active development.
