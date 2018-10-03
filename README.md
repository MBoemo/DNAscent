# Osiris
Software for detecting contiguous regions of base analogues incorporated in Oxford Nanopore reads.

## Downloading and Compiling Osiris
Clone the Osiris repository with the recursive flag so that the dependencies as well.
```shell
git clone --recursive https://github.com/MBoemo/Osiris.git
```
The Osiris directory will appear in your current directory.  Compile the software by running:
```shell
cd Osiris
make
```
This will put the Osiris executable into the Osiris/bin directory.

## Sample Workflow
We assume the following:
- you have a directory of 1D R9.5 450bp/s Oxford Nanopore fast5 reads (which may be in subdirectories) that you want to use for detection,
- these reads have been basecalled to fastq format using Albacore (version > 2.0),
- you have a reference/genome file for your reads in fasta format,
- you have aligned the reads fastq to a reference using an appropriate aligner and have an indexed BAM file.

Before running detection, the fast5 reads need to be indexed.  Run,
```shell
python Osiris/scripts/index.py -d path/to/fast5/reads
```
This will produce a file `index.osiris` which will be needed as input to the `Osiris detect` executable.

You can run `Osiris detect` (on 10 threads, for example) by running:
```shell
bin/Osiris detect -b /path/to/alignment.bam -r path/to/reference.fasta -m pore_models/BrdU_threshold2.0.model -i index.osiris -o output.detect -t 10
```
This will put a file called `output.detect` in your current directory.  Each read is given a header starting with a `>` character, followed by the read's ONT read ID, the chromosome this read aligned to, and the alignment's start and end on that chromosome.  The first column gives a position on the reference where a call was attempted, and the second column gives the log likelihood that the 6mer at this position contained a BrdU.  The third and fourth columns show the 6mer on the reference and the aligned 6mer on the basecalled read, respectively.  (Note that when high concentrations of analogues are present, this usually causes a significant disruption to Albacore's basecalling accuracy.)





## Dependencies
It relies on Penthus (https://github.com/MBoemo/Penthus) for its hidden Markov model functionality.  Development was done using gcc 5.4.0 on an Ubuntu 16.04 platform.
We have tried to keep dependencies to a bare minimum and only use standard software when dependencies were absolutely required.  The following are required for the scripts that prepare data for Osiris:
- Python 2.7 (https://www.python.org/downloads/),
- pysam (https://github.com/pysam-developers/pysam/releases),
- h5py (https://github.com/h5py/h5py/releases),
- samtools (https://github.com/samtools/samtools/releases),
- GraphMap (https://github.com/isovic/graphmap/releases).


## Tutorial Using Example Data
We have provided example data to be used in the following tutorial.  [Complete this later]

## Using C++ Osiris For Your Own Project
We assume the following:

- (for training only) you sequenced hairpin primers that have {NBNNNNN, NNNBNNN, NNNNNBN} domains along with their reverse complement domains, possibly pooled (and suitably barcoded), 
- you have a directory of 1D R9.5 450bp/s Oxford Nanopore reads (which may be in subdirectories) that you want to use for training or detection,
- these reads have been basecalled to fast5 format using Albacore (version > 2.0).

### Prepping Your Data
Before using Osiris, you have to prepare your data into a form that Osiris understands.  This step makes training easy to repeat, as you'll only have to open the fast5 files once.  This is done with a collection of provided Python scripts which are located in the scripts subdirectory.  All of these scripts can be run with:

```python
python <script>.py -h
```
to see the required and optional arguments.

demultiplexHairpin.py (for training) and demultiplex.py (for detection) do the following:
- goes through the reads directory, opens each fast5 file, and takes the basecalled sequence out of it,
- creates a reads.fasta file, where each basecalled sequence is named by the full path of to the fast5 file it came from,
- aligns reads.fasta against a reference with GraphMap,
- splits the alignments.sorted.bam file into separate .bam files for each reference.  

prepHairpinTrainingData.py does the following:
- takes one of the bam files produced by demultiplex.py and a reference fasta file that contains only the reference sequence for that bam file,
- groups the reads into bins, where each bin an identical NNNBNNN (or NBNNNNN, or NNNNNBN) domain,
- for bins that have the minimum number of training reads specified (using the -n flag), opens each fast5 file and converts the raw signal to picoAmps,
- creates a .foh file that is used as input for the Osiris executable.

### Training

Model training is done with the Osiris train executable.  It has the following required arguments:
- -r,--reference: a path to a reference file in fasta format.  This file should contain only one sequence (if it contains more than one, only the first one will be used),
- -p,--position: the position of the analogue in the domain, where NBNNNNN corresponds to 1and2, NNNBNNN corresponds to 3and4, and NNNNNBN corresponds to 5and6,
- -d,--trainingData: the .foh file made by prepHairpinData.py that you want to use as training data,
- -o,--output: the output path and name to a model file that Osiris train should produce.

In addition, there are the following optional arguments:
- -t,--threads: number of threads to run on (default is 1).  Note that Osiris train is helped significantly by multithreading.
- -c,--clipping: restricts training to the region on the reference where the analogue domain is located.  Increaes speed by several orders of magnitude, but at a slight cost to accuracy.  Note that this option is not available if the analogue domain is very close (within 20 bases) to the 3' or 5' end of the reference.
- -l,--log-file: path to a log file that shows the training at each position on the reference for each bin (default is none).
 
An example command is:
```shell
Osiris train -r path/to/reference.fasta -p 3and4 -d path/to/trainingData.foh -o path/to/outputModelFile.model -t 20 -c -l path/to/logFile.log
```
This command takes a reference file called reference.fasta and training data from trainingData.foh, trains for a base analogue at positions 3 and 4 in a 6mer, and writes the results to outputModelFile.model and logFile.log.  It does this by multithreading on 20 threads and using clipping.

As previously mentioned, the log file (should you choose to specify one) will show the training at each position on the reference for each bin.  A snippet of the log file might look like this:
```c++
>CCATTCG
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
The header following the > character shows the 7mer for this bin (the values of N in NBNNNNN, NNNBNNN, or NNNNNBN).  The columns, from left to right, are:
- position on the reference,
- the 6mer at that position,
- original mean from the Oxford Nanopore model file (before training),
- trained mean (after training),
- original standard deviation from the Oxford Nanopore model file (before training),
- trained standard deviation (after training).

### Detection

[Complete this later]
