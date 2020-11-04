//----------------------------------------------------------
// Copyright 2019-2020 University of Oxford
// Written by Michael A. Boemo (mb915@cam.ac.uk)
// This software is licensed under GPL-3.0.  You should have
// received a copy of the license with this software.  If
// not, please Email the author.
//----------------------------------------------------------

#ifndef SENSE_H
#define SENSE_H

#include <cassert>
#include <vector>
#include "error_handling.h"

/*function prototypes */
int sense_main( int argc, char** argv );

class DetectedRead{

	public:
		std::vector< unsigned int > positions;
		std::vector< double > brduCalls;
		std::string readID, chromosome, strand, header;
		int mappingLower, mappingUpper;
		std::vector<std::vector<float>> probabilities;
		std::vector<std::pair<int,int>> origins;
		std::vector<std::pair<int,int>> terminations;
		std::vector<float> tensorInput;
		void trim(unsigned int trimFactor){
			
			assert(positions.size() > trimFactor and brduCalls.size() > trimFactor and positions.size() == brduCalls.size());
			unsigned int cropFromEnd = positions.size() % trimFactor;
			brduCalls.erase(brduCalls.end() - cropFromEnd, brduCalls.end());
			positions.erase(positions.end() - cropFromEnd, positions.end());
			assert(positions.size() % trimFactor == 0 and brduCalls.size() % trimFactor == 0);
		}
		void generateInput(void){

			for (size_t i = 0; i < positions.size(); i++){
				tensorInput.push_back(brduCalls[i]);
			}
		}
};

#endif

