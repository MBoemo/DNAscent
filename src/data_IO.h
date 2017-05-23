//----------------------------------------------------------
// Copyright 2017 University of Oxford
// Written by Michael A. Boemo (michael.boemo@path.ox.ac.uk)
//----------------------------------------------------------


#include <fstream>
#include <sstream>
#include <string>
#include <iostream>
#include <algorithm>
#include <map>
#include "data_structures.h"


/*function prototypes */
std::string import_reference( std::string );
std::map< std::string, std::pair< double, double > > import_poreModel( std::string );
std::map< std::string, std::vector< std::vector< double > > > import_foh( std::string );
std::vector< detectionTuple > import_fdh( std::string & );
