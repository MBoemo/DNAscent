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


struct detectionTuple{

	std::string filename;
	std::string basecalls;
	std::vector< double > events;

};
