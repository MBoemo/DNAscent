# Osiris
Python/C++ HMM package for detecting base analogues in Oxford Nanopore reads.  It works by creating a branching HMM from a reference file.  One branch is for 6-mers that only contain thymidine, and the other branch is for 6-mers that contain at least one analogue.  The Python flavour of Osiris uses some of the HMM libraries from pomegranate (https://github.com/jmschrei/pomegranate).  Instructions for using Osiris in Python are as follows.

## Python Flavour
First install pomegranate (https://github.com/jmschrei/pomegranate), which is a HMM library written implemented in Cython that serves as a backend for Osiris.

### Installation and Setup
Install Osiris by running:
```python
python setup.py install
```      
Import Osiris functions with:
```python
import Osiris as Osiris
```

### Basic Usage
The first step is to build a HMM from a reference file.  If the DNA samples include base analogues, then this HMM will have a branching structure.  The inputs are:
(1.) a reference .fasta file
(2.) a base model, which has the emission data (mean, standard deviation) for all 6-mers, formatted for Nanopolish,
(3.) (optional) an Osiris BaseAnalogue object that contains the analogue's name as well as data for the concentration and emission probabilities.

The first argument is trivial, and the second argument can be found by using the same 6-mer model that Nanopolish uses.  To get the third argument, simply create the analogue object.  

```python
analogue = Osiris.BaseAnalogue()
analogue.set_concentration(0.50)
analogue.set_emissions(path-to-nanopolish-eventalign-file)
```
Note that the position of the analogue is indexed starting from 1 (unlike normal Python indexing which starts from 0).

If we're not interested in analogues, you can build the HMM by running:
```python
model = Osiris.build_hmm('reference.fasta','modelfile.model')
```
To build a HMM with branching structure for a base analogue, run:
```python
model = Osiris.build_hmm('reference.fasta',modelfile.model',analogue)
```
In the HMM, the reference file determines the branching structure where the HMM is in the branched hidden states if there is an analogue somewhere in the 6-mer.    

NOTE: The HMM is optimised and error checked using pomegranate's built-in optimisation routine, which can be very slow.  These optimisation routines ensure that the leaving probabilities for all states sum to 1, there are no needless silent states with a single exit transition of probability 1, and that there are no "orphan" states.  By default, Osiris does not use these optimisation/checking routines, but instead has been subjected to a large number of tests to ensure it builds the model correctly.  If you make changes to the source code, it is strongly recommended that you re-enable the optimisation/checking.

### Model Training
Import sequences from Nanopolish eventalign files with:
```python
trainingData = Osiris.import_sequences('training_data.eventalign')
```
You can run Baum-Welch training on the model with pomegranate's fit function.  This has options for parallel processing and conversion thresholds.  You can perform 30 training iterations on 20 cores by running:
```python
model.fit(trainingData,max_iterations=30,n_jobs=20)
```

### Saving/Loading Model Objects
Building and training the model are both computationally expensive, but the actual model objects tend to be quite small.  To prevent unnecessary re-building or re-training in the result of a crash or broken pipe, Osiris automatically exports built and trained models as json files to the working directory.  In general, these files take up negligible disk space (less than 15 MB).  To reload these saved files, do the folling:
```python
import json
import pomegranate as pm
f = open('model.json','r')
loaded_json = json.load(f)
loaded_model = pm.HiddenMarkovModel.from_json(loaded_json)
```
Writing and loading with json is on the order of seconds, whereas pickle and cPickle will be on the order of hours.  N.B. by default, pomegranate's from_json function will re-optimise the model with full optimisation.  It's worth turning that off.

##C++ Flavour
Under development.
