//----------------------------------------------------------
// Copyright 2017 University of Oxford
// Written by Michael A. Boemo (michael.boemo@path.ox.ac.uk)
// This software is licensed under GPL-3.0.  You should have
// received a copy of the license with this software.  If
// not, please Email the author.
//----------------------------------------------------------

#ifndef BUILD_MODEL_H
#define BUILD_MODEL_H

#include <fstream>
#include <sstream>
#include <string>
#include <iostream>
#include <algorithm>
#include <map>


/*function prototypes */
std::stringstream buildAndTrainHMM( std::string &, std::map< std::string, std::pair< double, double > > &, std::vector< std::vector< double > > &, int &, bool );
double buildAndDetectHMM( std::string &, std::map< std::string, std::pair< double, double > > &, std::map< std::string, std::pair< double, double > > &, std::vector< double > &, bool );

#endif
