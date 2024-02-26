//----------------------------------------------------------
// Copyright 2019-2020 University of Oxford
// This software is licensed under GPL-3.0.  You should have
// received a copy of the license with this software.  If
// not, please Email the author.
//----------------------------------------------------------

#ifndef EVENT_HANDLING_H
#define EVENT_HANDLING_H

#include "data_IO.h"

struct eventDataForRead {

	std::vector< event > events;
	std::vector< std::pair< unsigned int, unsigned int > > eventAlignment;
	std::map<unsigned int, double> posToScore;
	bool failed = false;
	double qualityScore;
};

void normaliseEvents( read &, bool );
void bulk_getEvents( std::string fast5Filename, std::string readID, std::vector<double> &raw );
void getEvents( std::string fast5Filename, std::vector<double> &raw );

#endif
