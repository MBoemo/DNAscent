# Osiris
Python/C++ software for detecting base analogues in Oxford Nanopore reads.  The two general uses are (1) determining the current across the pore that is produced by a 6mer that contains a single base analogue, and (2) determining where base analogues are incorporated in a nanopore read.  The Python flavour of Osiris uses some of the HMM libraries from pomegranate (https://github.com/jmschrei/pomegranate) while the C++ flavour of Osiris uses Penthus (https://github.com/MBoemo/Penthus).  Instructions for using Osiris in Python are as follows.

## Python Flavour
Dependencies:
- samtools (http://samtools.sourceforge.net/),
- pysam (https://github.com/pysam-developers/pysam),
- pomegranate (https://github.com/jmschrei/pomegranate).

Pomegranate is a Cython library of hidden markov model (HMM) algorithms that serves as a backend for Osiris.  All of these dependencies should be installed before installing Osiris.

### Installation and Setup
Install Osiris by running:
```python
python setup.py install
```      
Import Osiris functions with:
```python
import Osiris as osi
```

### Finding Base Analogue Emissions
The general inputs for Osiris are:
- A reference fasta file, which should contain the reference sequence of all samples used in the run.  They should each have unique names.  The only admissable characters in the reference sequences are A, T, G, C, and N.  Note: if there is a base analogue, it should be replaced by the standard base in the reference.  For example, BrdU (a thymidine analogue) should be replaced by a T in the reference file.
- A directory that contains all of the fast5 files from the run.  The fast5 files can be in subdirectories of this directory.

A nanopore run may contain several samples that were barcoded using the ONT barcoding kit.  The first step is to sort these reads by the reference that they best align to and perform some quality control on the reads.  This is done by the Osiris function alignAndSort.
```python
osi.alignAndSort('full-path-to-top-level-reads-directory','full-path-to-reference-fasta',number-of-threads)
```
Full paths to the reads and reference should be provided as strings, and the third argument specifies the number of threads on which to run BWA-MEM alignment.  The function will export all the fast5 reads to a fasta file, perform a sequence alignment against the reference, ignore reads that are reverse complements or have below 80% coverage of one of the references, and finally sort the BAM file into a separate BAM file for each reference.  The current working directory will now contain reads.fasta and a BAM file for each reference.

Once the function finishes, create a new reference fasta file that contains only the reference of interest (the one that contains an analogue).  Import the new reference file by running:

```python
reference = osi.import_reference('mySingleEntryReference.fasta')
```

Osiris accepts two types of training data.  The first has a base analogue (BrdU, for example, which we will denote by B) in a fixed context, such as 5'-ATGCCBTCCAG-3'.  The second has an analogue flanked by three N's with its reverse complement downstream, such as 5'-ATCG...NNNBNNN...NNNANNN...-3'.  We shall follow a workflow for the latter type of training data.  The case for the former is similar.

An example of a command that imports training data is:
```python
trainingData = osi.import_HairpinTrainingData(reference,'mySingleEntryReference.bam','template_median68pA.5mers.model',location-of-redundant-A,minimum-reads-threshold)
```
Here, the first input is the output of the previous step.  The alignment file mySingleEntryReference.bam is the BAM file created by the alignAndSort function: it is the BAM file that corresponds to the reference.  The 5mer model file can be found in the pore_models Osiris directory, and is required here to normalise for the Metrichor shift and scale parameters.  You can either copy it to your current working directory, or put the full path as the argument.  The input location-of-redundant-A is an integer that gives the location in the reference of the A base of the NNNANNN domain.  Please note that in Osiris, locations are indexed from zero (as is typical in Python) so that the first base in the reference has index 0.  Finally, the last argument is an integer that gives a cutoff for the minimum number of reads per 6mer that Osiris will use as training data.  So if we enter 20, and your run only has less than 20 reads for a 6mer, Osiris will not incorporate that 6mer into the new pore model because it is considered to have too little data.

To train the model, run:
```python
trainedEmissions = osi.trainForContextAnalogue(trainingData, reference, 'template_median68pA.model', threads)
```
This function builds a HMM for the reference, and uses the training data created in the previous step to adjust the emission parameters to those of the base analogue.  The inputs {trainingData, reference} were generated in previous steps.  The model file template_median68pA.model is provided in the pore_models directory.  Note that this is a 6mer model, as opposed to the 5mer model that was used to normalise the training data.  The function must be provided with the number of cores that model training is allowed to run on.

Finally, to visualise the new base analogue emission means and standard deviations, Osiris will build a new pore model file.  Run:
```python
osi.export_poreModel(trainedEmissions, 'output-pore-model-file.model')
```
The input trainedEmissions was created in the prevous step, and 'BrdU_emissions.model' is the filename of the pore model that you want the function to create.  After calling export_poreModel, you should see the pore model file appear in your working directory.

### Calling Base Analogues
The second main use of Osiris is using trained analogue emissions to detect base analogues that have been randomly incorporated in a Nanopore read.  Assume we have a trainedEmissions dictionary, which could have been made by following the steps above or by importing an Osiris base analogue pore model file.  Create a base analogue with the analogue emissions and the concentration by running:
```python
import Osiris as osi
BrdU = osi.BaseAnalogue(trainedEmissions,0.2)
```
Here, the base analogue (BrdU) has replaced 20% of all T's in the read.  Build an Osiris hidden markov model by running:
```python
hmm = osi.build_RandIncHMM('reference.fasta','template_median68pA.model',analogue=BrdU)
```
The 6mer template model file can be found in the pore_models directory, and the refernece fasta file should only include one reference sequence for the reads you want to analyse.

Suppose we want to check read.fast5 to determine where BrdU has been incorporated.  Normalise the events for shift and scale with,
```python
normalisedEvents = osi.calculate_normalisedEvents(['read.fast5'], 'template_median68pA.5mers.model')
```
Find the locations of BrdU in the read by running:
```python
BrdU_positions = osi.callAnalogue(hmm, normalisedEvents[0])
```
The variable BrdU_positions is a list of positions in the reference where Osiris has called a base analogue.

### Saving/Loading Objects
Normalising the training data and building reference-specific HMMs can be computationally expensive, but the space required by both training data and model objects tends to be quite small.  To prevent having to recompute these, you can save the training data or model to a file using JSON.  To save and load training data, run the following:
```python
import json

#export
f = open('trainingData.json','w')
json.dump(trainingData,f)
f.close()

#import
f = open('trainingData.json','r')
trainingData = json.load(f)
f.close()
```
To do the equivalent for a HMM object, run:
```
import json
import pomegranate as pm

#export
hmm_json = hmm.to_json()
f = open('model.json','w')
json.dump(hmm_json,f)
f.close()

#import
f = open('model.json','r')
loaded_json = json.load(f)
hmm = pm.HiddenMarkovModel.from_json(loaded_json)
```
Writing and loading with JSON is on the order of seconds, whereas pickle and cPickle will be on the order of hours.  Please note that by default, pomegranate's from_json function will re-optimise the model with full optimisation.  It's worth turning that off by editing the code.

## C++ Flavour
The C++ flavour of Osiris uses Penthus as its machine learning core.  You can compile Osiris by navigating to the Osiris directory and running:
```shell
make
```
This will build the main Osiris executable in the Osiris/bin directory.  To train a base analogue pore model on hairpin primer training data, run:
```shell
./Osiris train -r /full/path/to/reference.fasta -m /full/path/to/template_median68pA.model -d /full/path/to/trainingData.foh -o /full/path/to/desiredOutputFile.txt
```
You can also use the optional -t flag, which specifies the number of threads for parallel processing.  The output file should look like:
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

