#!/usr/bin/env python

#----------------------------------------------------------
# Copyright 2017 University of Oxford
# Written by Michael A. Boemo (michael.boemo@path.ox.ac.uk)
#----------------------------------------------------------

from distutils.core import setup
from distutils.extension import Extension
import numpy as np

filenames = [ 'build_model',
	      'data_IO',
	      'train',
	      'utility'
            ]

modules= ["Osiris/{}".format(name) for name in filenames]

setup(
      name='Osiris',
      version='1.0',
      description='Software for detecting BrdU and base analogues in Oxford Nanopore reads.',
      author='Michael A. Boemo',
      author_email='michael.boemo@path.ox.ac.uk',
      url='https://github.com/MBoemo/osiris',
      license='LICENSE',
      packages=['Osiris'],
      py_modules=modules
     )
