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
git checkout v1.0.0
make
```
This will put the DNAscent executable into the DNAscent/bin directory.  A typical compile time for DNAscent and its dependencies is 5 minutes.

## Documentation
Please see the [documentation](https://dnascent.readthedocs.io) for detailed usage instructions, descriptions of DNAscent's subprograms, and an example workflow.

## Citation
Muller CA, Boemo MA, Spingardi P, Kessler BM, Kriaucionis S, Simpson JT, Nieduszynski CA. Capturing the dynamics of genome replication on individual ultra-long nanopore sequence reads. *Nature Methods* 2019;16:429-436. [[Journal Link](https://www.nature.com/articles/s41592-019-0394-y)]

## Questions and Bugs
Should any bugs arise or if you have basic usage questions, please raise a [GitHub issue](https://github.com/MBoemo/DNAscent/issues). For more detailed discussions or collaborations, please Email Michael Boemo at mb915@cam.ac.uk.
