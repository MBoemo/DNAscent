#Osiris
Software for detecting base analogues in Oxford Nanopore reads.  The two general uses are: (1.) determining the characteristic current produced by a 6mer that contains a single base analogue, and (2.) determining where base analogues are incorporated in individual Oxford Nanopore reads.  The Python flavour of Osiris uses some of the HMM libraries from pomegranate (https://github.com/jmschrei/pomegranate) while the C++ flavour of Osiris uses Penthus (https://github.com/MBoemo/Penthus).

##Dependencies
We have tried to keep dependencies to a bare minimum and only use standard software when dependencies were absolutely required.  The following are required for the scripts that prepare data for Osiris:
- Python 2.7 (https://www.python.org/downloads/),
- pysam (https://github.com/pysam-developers/pysam/releases),
- h5py (https://github.com/h5py/h5py/releases),
- NumPy (https://github.com/numpy/numpy/releases),
- samtools (https://github.com/samtools/samtools/releases),
- GraphMap (https://github.com/isovic/graphmap/releases).

## Getting and Compiling Osiris
Clone the Osiris repository with the recursive flag so that you get Penthus as well.
```shell
git clone --recursive https://github.com/MBoemo/Osiris.git
```
You can compile C++ Osiris by running:
```shell
make
```
This will put the Osiris executable into Osiris/bin.  Optionally, you can add the Osiris executable to your path with:
```shell
PATH=$PATH:/path/to/bin
export PATH
```

##Tutorial Using Example Data
We have provided example data to be used in the following tutorial.  [Complete this later]

##Using C++ Osiris For Your Own Project
We assume the following:

- (for training only) you sequenced hairpin primers that have {NBNNNNN, NNNBNNN, NNNNNBN} domains along with their reverse complement domains, 
- you have a directory of 1D R9.5 450bp/s Oxford Nanopore reads (which may be in subdirectories) that you want to use for training or detection,
- these reads have been basecalled to fast5 format using Albacore >v1.2 with event detection.

###Prepping Your Data
Before using Osiris, you have to prepare your data into a form that Osiris understands.  This is done with a collection of provided Python scripts which are located in the scripts subdirectory.  All of these scripts can be run with

```python
python <script>.py -h
```
to see the required and optional arguments.

demultiplex.py (for detection) and demultiplexHairpin.py (for training) look through the reads directory, create a reads.fasta file, performs an alignment using GraphMap, and then splits the alignments.sorted.bam file into separate .bam files for each reference.  

prepHairpinTrainingData.py takes one of the .bam files from demultiplex.py and a reference .fasta file that contains only the reference sequence for the .bam file and outputs a .foh file of normalised events.  This .foh file is used as the input for the Osiris executable.  Note that prepHairpinTrainingData.py requires a 5mer model.  This 5mer model file is provided in the pore_models subdirectory.

###Training

To train a base analogue pore model on hairpin primer training data, run:
```shell
./Osiris train -r /full/path/to/reference.fasta -p 3and4 -om /full/path/to/pore_models/template_median68pA.6mer.model -d /full/path/to/trainingData.foh -o /full/path/to/outputModelFile.trained -t 20 -sc 35
```
You can specify a log file with the -l flag.  In addition to the output model file, the log file will give the trained emission at each position.  The log file should look like:
```c++
>CCAATCG
state	info	oriMu	trMu	oriSig	trSig
0_M1	AATGTA	107.047	101.258	3.28075	6.2622
1_M1	ATGTAC	80.9837	72.0795	2.32535	2.39806
2_M1	TGTACT	97.6125	97.8708	2.27356	3.54221
3_M1	GTACTT	85.0232	85.4749	1.48386	0.89079
4_M1	TACTTC	81.5373	80.7603	1.57174	1.13429
5_M1	ACTTCG	95.3118	93.7051	1.66458	2.20321
6_M1	CTTCGT	103.454	102.039	2.38881	1.34545
7_M1	TTCGTT	89.9279	90.747	2.01772	2.10616
8_M1	TCGTTC	62.4801	62.7164	1.84319	2.16708
9_M1	CGTTCA	93.5076	91.8787	2.18633	1.82847
```
The columns, from left to right, are:

- position on the reference,
- 6mer for that state,
- original mean (before training),
- trained mean (after training),
- original standard deviation (before training),
- trained standard deviation (after training).

###Detection

[Complete this later]

## Python Flavour
NOTE: The Python flavour of Osiris is deprecated; we strongly recommend that you use the C++ flavour.  

Pomegranate is a Cython library of hidden Markov model (HMM) algorithms that serves as a backend for Python Osiris.  Clone and install it from (https://github.com/jmschrei/pomegranate), then clone the Osiris repository.  After the repository has been cloned into your working directory, navigate to the Osiris directory.  Install Python Osiris by running:
```python
python setup.py install
```
Note that you may have to install Python Osiris as an administrator, depending on your system setup.  Import the Osiris functions from a Python script or shell with:
```python
import Osiris as osi
osi.some_Osiris_function
```
