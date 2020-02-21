//----------------------------------------------------------
// Copyright 2019 University of Oxford
// Written by Michael A. Boemo (mb915@cam.ac.uk)
// This software is licensed under GPL-2.0.  You should have
// received a copy of the license with this software.  If
// not, please Email the author.
//----------------------------------------------------------

#ifndef DETECT_H
#define DETECT_H

#include <utility>
#include <string>
#include <vector>
#include "data_IO.h"
#include "../htslib/htslib/hts.h"
#include "../htslib/htslib/sam.h"

/*function prototypes */
int detect_main( int argc, char** argv );

std::vector< unsigned int > getPOIs( std::string &, int );
void countRecords( htsFile *, hts_idx_t *, bam_hdr_t *, int &, int, int );
void parseIndex( std::string, std::map< std::string, std::string > &, bool & );
void parseCigar(bam1_t *, std::map< unsigned int, unsigned int > &, int &, int & );
std::string getQuerySequence( bam1_t * ); 
double sequenceProbability( std::vector <double> &,
				std::string &,
				size_t,
				bool,
				PoreParameters,
				int );

#endif
