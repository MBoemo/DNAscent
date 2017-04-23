//----------------------------------------------------------
// Copyright 2017 University of Oxford
// Written by Michael A. Boemo (michael.boemo@path.ox.ac.uk)
//----------------------------------------------------------


#include "data_IO.h"
#include "../Penthus/src/hmm.h"


/*function prototypes */
std::stringstream buildAndTrainHMM( std::string &, std::map< std::string, std::pair< double, double > > &, std::vector< std::vector< double > > & );
