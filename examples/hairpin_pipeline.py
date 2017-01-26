#!/usr/bin/env python


import Osiris as osi
import json

#import reference file
reference = osi.import_reference('reference.fasta')

#import and format training data
trainingData = osi.import_HairpinTrainingData('pooled_BrdU_reads.fasta', reference, 'merged.bam', 'template_median68pA.5mers.model', 113, 20)

#export the training data as json so we can import it later if we want
f = open('trainingData_BrdU_merged.json','w')
json.dump(trainingData,f)
f.close()

#train the HMM on the analogue data
trainedEmissions = osi.trainForContextAnalogue(trainingData, reference, 'template_median68pA.model', 0.5, 20)

#create a pore model file to see the mean and standard deviation of the analogue 6mers' signal
osi.export_poreModel(trainedEmissions, 'BrdU_emissions.model')
