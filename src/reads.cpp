//----------------------------------------------------------
// Copyright 2019 University of Oxford
// This software is licensed under GPL-3.0.  You should have
// received a copy of the license with this software.  If
// not, please Email the author.
//----------------------------------------------------------

#include <string>
#include <map>
#include <vector>
#include <iostream>
#include <cassert>
#include <memory>
#include "reads.h"

std::vector< std::vector<DNAscent::read *> > sortReadsByFilename(std::vector<DNAscent::read> &buffer){

	std::map<std::string, std::vector<DNAscent::read *> > filenameToReads;
	for (size_t i = 0; i < buffer.size(); i++){
	
		if (filenameToReads.count(buffer[i].filename) == 0){
		
			filenameToReads[buffer[i].filename] = {&buffer[i]};
		}
		else{
		
			filenameToReads[buffer[i].filename].push_back(&buffer[i]);
		}
	}
	
	std::vector< std::vector<DNAscent::read *> > batchedReads;
	for (auto i = filenameToReads.begin(); i != filenameToReads.end(); i++){
	
		batchedReads.push_back( i -> second );
	}
	
	return batchedReads;
}
