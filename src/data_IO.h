//----------------------------------------------------------
// Copyright 2017 University of Oxford
// Written by Michael A. Boemo (michael.boemo@path.ox.ac.uk)
//----------------------------------------------------------


#include <fstream>
#include <string>
#include <iostream>
#include <algorithm>
#include <map>


/*function prototypes */
std::string import_reference( std::string );
std::map< std::string, std::pair< double, double > > import_poreModel( std::string );
