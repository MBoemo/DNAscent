//----------------------------------------------------------
// Copyright 2017 University of Oxford
// Written by Michael A. Boemo (michael.boemo@path.ox.ac.uk)
// This software is licensed under GPL-2.0.  You should have
// received a copy of the license with this software.  If
// not, please Email the author.
//----------------------------------------------------------

#ifndef PSL_H
#define PSL_H

#include <vector>

struct Track{

	int lowerBound, upperBound;
};

struct readDetection{

	std::vector< int > positions;
	std::vector< double > BrdUProb;
	std::string readID, chromosome;
	int mappingLower, mappingUpper;
	std::vector< Track > tracks;
	std::string direction;
};


/*function prototypes */
int psl_main( int argc, char** argv );

#endif
