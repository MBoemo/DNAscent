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


/*function prototypes */
std::stringstream buildAndTrainHMM( std::string &, std::map< std::string, std::pair< double, double > > &, std::vector< std::vector< double > > &, int & );
double buildAndDetectHMM( std::string &, std::map< std::string, std::pair< double, double > > &, std::map< std::string, std::pair< double, double > > &, std::vector< double > &, int &, bool );
