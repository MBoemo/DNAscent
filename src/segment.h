//----------------------------------------------------------
// Copyright 2020 University of Cambridge
// Written by Michael A. Boemo (mb915@cam.ac.uk)
// This software is licensed under GPL-3.0.  You should have
// received a copy of the license with this software.  If
// not, please Email the author.
//----------------------------------------------------------

#ifndef SEGMENT_H
#define SEGMENT_H

#include <cassert>
#include <vector>
#include <iostream>
#include <fstream>
#include "error_handling.h"

int segment_main( int argc, char** argv );

struct segmentArgs {

	std::string detectFilename;
	int minLength = 1000;
	unsigned int threads = 1;
};


#endif

