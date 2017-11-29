//----------------------------------------------------------
// Copyright 2017 University of Oxford
// Written by Michael A. Boemo (michael.boemo@path.ox.ac.uk)
// This software is licensed under GPL-3.0.  You should have
// received a copy of the license with this software.  If
// not, please Email the author.
//----------------------------------------------------------


#ifndef COMMON_H
#define COMMON_H


#include <algorithm>
#include <vector>
#include <utility>
#include <iterator>
#include <map>
#include <sstream>
#include <iostream>


inline std::string reverseComplement( std::string DNAseq ){

	std::reverse( DNAseq.begin(), DNAseq.end() );
	std::string revComp;

	for ( std::string::iterator i = DNAseq.begin(); i < DNAseq.end(); i++ ){
	
		switch( *i ){
			case 'A' :
				revComp += 'T';
				break;
			case 'T' :
				revComp += 'A';
				break;
			case 'G' :
				revComp += 'C';
				break;
			case 'C' :
				revComp += 'G';
				break;
			default:
				std::cout << "Exiting with error.  Invalid character passed to reverse complement function.  Must be A, T, G, or C." << std::endl;
				exit( EXIT_FAILURE );
		}
	}
	return revComp;
}


/*function prototypes */
void displayProgress( int, int );
std::vector< std::string > split( std::string, char );
int argMin( std::vector< double > & );


#endif
