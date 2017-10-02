//----------------------------------------------------------
// Copyright 2017 University of Oxford
// Written by Michael A. Boemo (michael.boemo@path.ox.ac.uk)
//----------------------------------------------------------

#ifndef DATA_IO_H
#define DATA_IO_H

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
void export_poreModel( std::map< std::string, std::vector< double > > &, std::string &);

#endif
