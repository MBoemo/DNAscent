# DNAscent
Software for detecting contiguous regions of base analogues incorporated in Oxford Nanopore reads.  Development was done using gcc 5.4.0 on an Ubuntu 16.04 platform.

## Downloading and Compiling DNAscent
Clone the DNAscent repository with the recursive flag so that the dependencies are cloned as well.
```shell
git clone --recursive https://github.com/MBoemo/DNAscent.git
```
The DNAscent directory will appear in your current directory.  Compile the software by running:
```shell
cd DNAscent
make
```
This will put the DNAscent executable into the DNAscent/bin directory.  A typical compile time for DNAscent and all of its dependencies is 5-7 minutes.

## Sample Workflow
We assume the following:
- you have a directory of 1D R9.5 450bp/s Oxford Nanopore fast5 reads (which may be in subdirectories) that you want to use for detection,
- these reads have been basecalled to fastq format using Albacore or Guppy (available from Oxford Nanopore),
- you have a reference/genome file for your reads in fasta format,
- you have aligned the reads fastq to a reference using an appropriate aligner and have an indexed BAM file.

Before running detection, the fast5 reads need to be indexed.  Run,
```shell
DNAscent index -d path/to/fast5/reads -s path/to/sequencing_summary.txt
```
(Note that both Albacore and Guppy produce a sequencing summary file.  While this is not a required argument, it will make indexing run much faster and is strongly recommended.)  This will produce a file `index.dnascent` which will be needed as input to the `DNAscent detect` executable.

You can run `DNAscent detect` (on 10 threads, for example) by running:
```shell
bin/DNAscent detect -b /path/to/alignment.bam -r path/to/reference.fasta -i index.DNAscent -o output.detect -t 10
```
This will put a file called `output.detect` in your current directory.  Each read is given a header starting with a `>` character, followed by the read's ONT read ID, the chromosome this read aligned to, and the alignment's start and end on that chromosome.  The first column gives a position on the reference where a call was attempted, and the second column gives the log likelihood that the 6mer at this position contained a BrdU.  The third and fourth columns show the 6mer on the reference and the aligned 6mer on the basecalled read, respectively.  (Note that when high concentrations of analogues are present, this usually causes a significant disruption to Albacore/Guppy's basecalling accuracy.)

At this point, you may wish to view the individual analogue calls genome-wide.  DNAscent can produce a PSL file on the results from `DNAscent detect` which can then be visualised in IGV or the UCSC Genome Browser.  Run
```shell
bin/DNAscent psl -d output.detect -r path/to/reference.fasta -o output
```
which will put a file `output.psl` in your working directory.

If you can approximate the probability at which cannonical bases are substituted for an analogue in an analogue-rich region, DNAscent can use the output from `DNAscent detect` to call regions of analogue incorporation.  Run,
```shell
bin/DNAscent regions -d output.detect -p 0.2 -o output.regions
```
where `p` is the probability that there is an analogue in any 6mer.  The file `output.regions` will have the same header for each read as in `output.detect`.  The first column specifies where the region started, the second column specifies where the region ended, and the third column gives a z-score specifying how well this region fit the model for analogue incorporation.  Scores near or above 0 indicate that this region fit the model well and is most likely a region of analogue incorporation.  Negative scores indicate that there were fewer analogue calls than would be expected in an analogue region.  The DNAscent regions executable can also call fork direction and origin location.  To enable this functionality, add the `--replication` flag when you call DNAscent regions.  This will produce a file `calledOrigins.dnascent` that lists the chromosomal coordinates of called origins on individual reads.

## Runtime
The runtime of DNAscent is linear with respect to the amount of data (in bases) passed to it.

## Dependencies
Cloning the repository recursively (see above) will provide all the required dependencies, so you don't need to find them yourself.  These are:
- pfasta (https://github.com/kloetzl/pfasta)
- fast5 (https://github.com/mateidavid/fast5.git)
- htslib (https://github.com/samtools/htslib.git)
- hdf5lib (https://support.hdfgroup.org/HDF5/)
- tinydir (https://github.com/cxong/tinydir.git)

Note that the high throughput sequencing library (htslib) requires bzlib and lzma for compression.  If you don't have these, apt-get lzma-dev and liblzma-dev.  In addition, pfasta requires libbsd on Linux.
