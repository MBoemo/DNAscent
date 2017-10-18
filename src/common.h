//----------------------------------------------------------
// Copyright 2017 University of Oxford
// Written by Michael A. Boemo (michael.boemo@path.ox.ac.uk)
// This software is licensed under GPL-3.0.  You should have
// received a copy of the license with this software.  If
// not, please Email the author.
//----------------------------------------------------------


#ifndef COMMON_H
#define COMMON_H


#include "math.h"
#include <algorithm>
#include <vector>
#include <utility>
#include <iostream>
#include <map>
#include <iterator>
#include <limits>
#include <sstream>


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
std::vector< std::vector< double > > filterEvents( std::string &, std::map< std::string, std::pair< double, double > > &, std::vector< std::vector< double > > & );
std::pair< int, int > subsequenceDynamicTimewarping( std::vector< double > &, std::vector< double > & );
std::vector< double > generateSignal( std::string &, std::map< std::string, std::pair< double, double > > & );
std::map< int, std::vector< int > > dynamicTimewarping( std::vector< double > &, std::vector< double > & );

#endif
