# Osiris
Python/C++ software for detecting base analogues in Oxford Nanopore reads.  The two general uses are (1) determining the current across the pore that is produced by a 6mer that contains a base analogue, and (2) determining where base analogues are incorporated in a nanopore read.  The Python flavour of Osiris uses some of the HMM libraries from pomegranate (https://github.com/jmschrei/pomegranate).  Instructions for using Osiris in Python are as follows.

## Python Flavour
Dependencies:
- samtools (http://samtools.sourceforge.net/),
- pysam (https://github.com/pysam-developers/pysam),
- pomegranate (https://github.com/jmschrei/pomegranate).

Pomegranate is a Cython library of hidden markov model (HMM) algorithms that serves as a backend for Osiris.

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
The inputs to Osiris are:
- A reference fasta file, which should contain the reference sequence of all samples used in the run.  They should each have unique names.  The only admissable characters in the reference sequences are A, T, G, C, and N.
- A directory that contains all of the fast5 files from the run.

The first step is to sort the reads by the reference that they align to and perform a quality control on the coverage of these reads.  The shell script prepData.sh is provided in the scripts directory for this purpose.  Paths to the reference file and reads directory should be put under the "inputs" header at the top of the script.  Save the script after editing, and run it with:
```shell
bash prepData.sh
```
The script will export the fast5 reads to a fasta file, perform a sequence alignment against the reference, ignore reads that have below 80% coverage of one of the references, and finally sorts the BAM file into a BAM file for each reference.

Once the script finishes, create a few reference fasta file that contains only the reference of interest (that contains an analogue).  Import the new reference file by running:

```python
reference = osi.import_reference('mySingleEntryReference.fasta')
```

Osiris accepts two types of training data.  The first has a base analogue (BrdU, for example, which we will denote by B) in a fixed context, such as 5'-ATGCCBTCCAG-3'.  The second has an analogue flanked by three N's with its reverse complement downstream, such as 5'-ATCG...NNNBNNN...NNNANNN...-3'.  We shall follow a workflow for the latter type of training data.  The case for the former is similar.

To import the training data, run:
```python
trainingData = import_HairpinTrainingData('reads.fasta','mySingleEntryReference.fasta','mySingleEntryReference.bam','template_median68pA.5mers.model',113,20)
```
Here, reads.fasta is the fasta file created by prepData.sh.  mySingleEntryReference.fasta is the same file that was used to build the reference, and mySingleEntryReference.bam is the BAM file created by prepData.sh that corresponds to the reference.  The 5mer model file can be found in the pore_models directory, and is required here to normalise for the Metrichor shift and scale parameters.  The input 113 is this example is the location on the reference of the A base in NNNANNN.  Please note that in Osiris, locations are indexed from zero, so the first base in the reference has index 0.  Finally, the last argument is the minimum number of reads that we're allowed to train on (20 in this case).

To train the model, run:
```python
trainedEmissions = osi.trainForContextAnalogue(trainingData, reference, 'template_median68pA.model', 0.5, 20)
```
This builds a HMM for the reference, and uses the training data to adjust the emission parameters to those of the base analogue.  The inputs trainingData and reference were generated in previous steps.  The model file template_median68pA.model is provided in the pore_models directory.  Note that this is a 6mer model, as opposed to the 5mer model that was used to normalise the training data.  The function must be provided with a convergence for Baum-Welch iterations (0.5 in this case) and the number of cores that model training is allowed to run on (20 in this case).

Finally, to visualise the new base analogue emission means and standard deviations, Osiris will build a new pore model file.  Run:
```python
osi.export_poreModel(trainedEmissions, 'BrdU_emissions.model')
```
The input trainedEmissions was created in the prevous step, and 'BrdU_emissions.model' is the filename of the pore model that you want the function to create.

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
'''
Find the locations of BrdU in the read by running:
```python
BrdU_positions = osi.callAnalogue(hmm, normalisedEvents[0])
```
The variable BrdU_positions is a list of positions in the reference where Osiris has called a base analogue.

### Saving/Loading Model Objects
Normalising the training data and training models can be computationally expensive, but the actual training data and model objects tend to be quite small.  To prevent having to recompute these, you can save the training data or model to a file using json.  To save and load training data, run the following:
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
Writing and loading with json is on the order of seconds, whereas pickle and cPickle will be on the order of hours.  Please note that by default, pomegranate's from_json function will re-optimise the model with full optimisation.  It's worth turning that off by editing the code.

##C++ Flavour
Under development.
