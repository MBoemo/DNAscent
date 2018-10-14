# DNAscent
Software for detecting contiguous regions of base analogues incorporated in Oxford Nanopore reads.  Development was done using gcc 5.4.0 on an Ubuntu 16.04 platform.

## Downloading and Compiling DNAscent
Clone the DNAscent repository with the recursive flag so that the dependencies as well.
```shell
git clone --recursive https://github.com/MBoemo/DNAscent.git
```
The DNAscent directory will appear in your current directory.  Compile the software by running:
```shell
cd DNAscent
make
```
This will put the DNAscent executable into the DNAscent/bin directory.

## Sample Workflow
We assume the following:
- you have a directory of 1D R9.5 450bp/s Oxford Nanopore fast5 reads (which may be in subdirectories) that you want to use for detection,
- these reads have been basecalled to fastq format using Albacore (version > 2.0),
- you have a reference/genome file for your reads in fasta format,
- you have aligned the reads fastq to a reference using an appropriate aligner and have an indexed BAM file.

Before running detection, the fast5 reads need to be indexed.  Run,
```shell
python DNAscent/scripts/index.py -d path/to/fast5/reads
```
(Note that this script requires h5py - see below.)  This will produce a file `index.dnascent` which will be needed as input to the `DNAscent detect` executable.

You can run `DNAscent detect` (on 10 threads, for example) by running:
```shell
bin/DNAscent detect -b /path/to/alignment.bam -r path/to/reference.fasta -m pore_models/BrdU_threshold2.0.model -i index.DNAscent -o output.detect -t 10
```
This will put a file called `output.detect` in your current directory.  Each read is given a header starting with a `>` character, followed by the read's ONT read ID, the chromosome this read aligned to, and the alignment's start and end on that chromosome.  The first column gives a position on the reference where a call was attempted, and the second column gives the log likelihood that the 6mer at this position contained a BrdU.  The third and fourth columns show the 6mer on the reference and the aligned 6mer on the basecalled read, respectively.  (Note that when high concentrations of analogues are present, this usually causes a significant disruption to Albacore's basecalling accuracy.)

At this point, you may wish to view the individual analogue calls genome-wide.  DNAscent can produce a PSL file on the results from `DNAscent detect` which can then be visualised in IGV or the UCSC Genome Browser.  Run
```shell
bin/DNAscent psl -d output.detect -r path/to/reference.fasta -o output
```
which will put a file `output.psl` in your working directory.

If you can approximate the probability at which cannonical bases are substituted for an analogue in an analogue-rich region, DNAscent can use the output from `DNAscent detect` to call regions of analogue incorporation.  Run,
```shell
bin/DNAscent regions -d output.detect -p 0.5 -o output.regions
```
where `p` is the probability that there is an analogue in any 6mer.  The file `output.regions` will have the same header for each read as in `output.detect`.  The first column specifies where the region started, the second column specifies where the region ended, and the third column gives a z-score specifying how well this region fit the model for analogue incorporation.  Scores near 0 indicate that this region fit the model well and is most likely a region of analogue incorporation.  Negative scores indicate that there were fewer analogue calls than would be expected in an analogue region.

## Dependencies
Cloning the repository recursively (see above) will provide all the required dependencies.  The list of dependencies is as follows:
- h5py(https://github.com/h5py/h5py.git)
- Penthus (https://github.com/MBoemo/Penthus.git)
- fast5 (https://github.com/mateidavid/fast5.git)
- htslib (https://github.com/samtools/htslib.git)
- hdf5lib (https://support.hdfgroup.org/HDF5/)
