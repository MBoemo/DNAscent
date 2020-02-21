//----------------------------------------------------------
// Copyright 2019 University of Oxford
// Written by Michael A. Boemo (mb915@cam.ac.uk)
// This software is licensed under GPL-2.0.  You should have
// received a copy of the license with this software.  If
// not, please Email the author.
//----------------------------------------------------------

#ifndef EVENT_HANDLING_H
#define EVENT_HANDLING_H

#include "data_IO.h"

struct eventDataForRead {

	std::vector< double > normalisedEvents;
	std::vector< std::pair< unsigned int, unsigned int > > eventAlignment;
	std::map<unsigned int, double> posToScore;
	bool failed = false;
	double qualityScore;
};

void normaliseEvents( read & );
void bulk_getEvents( std::string fast5Filename, std::string readID, std::vector<double> &raw );
void getEvents( std::string fast5Filename, std::vector<double> &raw );

#endif
