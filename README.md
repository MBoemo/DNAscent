# DNAscent
[![Documentation Status](https://readthedocs.org/projects/dnascent/badge/?version=latest)](https://dnascent.readthedocs.io/en/latest/?badge=latest)

Software for detecting contiguous regions of base analogues incorporated in Oxford Nanopore reads.  Development was done using gcc 5.4.0 on an Ubuntu 16.04 platform.

## Downloading and Compiling DNAscent
The recommended OS for DNAscent is Ubuntu 16.04.  Clone the DNAscent repository with the recursive flag so that the dependencies are cloned as well.
```shell
git clone --recursive https://github.com/MBoemo/DNAscent.git
```
The DNAscent directory will appear in your current directory.  Switch to the latest tagged version and compile the software by running:
```shell
cd DNAscent
git checkout 3.1.2
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
