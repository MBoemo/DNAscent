//----------------------------------------------------------
// Copyright 2019-2020 University of Oxford
// This software is licensed under GPL-3.0.  You should have
// received a copy of the license with this software.  If
// not, please Email the author.
//----------------------------------------------------------

#ifndef DATA_IO_H
#define DATA_IO_H


#include <string>
#include <map>
#include <vector>
#include <iostream>
#include <cassert>
#include <unordered_map>
#include <omp.h>


struct IndexEntry{
	size_t batchIndex;
	size_t rowIndex;
	std::string filepath;
};


std::map< std::string, std::string > import_reference( std::string );
std::map< std::string, std::string > import_reference_pfasta( std::string );
std::vector< std::pair< double, double > > import_poreModel_staticStdv( std::string, unsigned int);
std::vector< std::pair< double, double > > import_poreModel_fitStdv( std::string, unsigned int);
std::string getExePath(void);
std::string getGitCommit(void);
unsigned int kmer2index(std::string &, unsigned int);
void parseIndex( std::string, std::map< std::string, IndexEntry > & );

#endif
