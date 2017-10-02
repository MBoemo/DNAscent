//----------------------------------------------------------
// Copyright 2017 University of Oxford
// Written by Michael A. Boemo (michael.boemo@path.ox.ac.uk)
// This software is licensed under GPL-3.0.  You should have
// received a copy of the license with this software.  If
// not, please Email the author.
//----------------------------------------------------------

#ifndef DATA_STRUCTURES_H
#define DATA_STRUCTURES_H

#include <fstream>
#include <sstream>
#include <string>
#include <iostream>
#include <algorithm>
#include <map>


struct detectionTuple{

	std::string filename;
	std::string basecalls;
	std::vector< double > events;
};

#endif
