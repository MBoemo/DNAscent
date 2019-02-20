//----------------------------------------------------------
// Copyright 2019 University of Oxford
// Written by Michael A. Boemo (michael.boemo@path.ox.ac.uk)
// This software is licensed under GPL-2.0.  You should have
// received a copy of the license with this software.  If
// not, please Email the author.
//----------------------------------------------------------


#include <random>
#include <limits>
#include "math.h"
#include "common.h"
#include "error_handling.h"


void displayProgress( int current, int total ){
/*poached from https://stackoverflow.com/questions/14539867/how-to-display-a-progress-indicator-in-pure-c-c-cout-printf */

	double progress = (double) current / (double) total;
	int barWidth = 70;

	if ( progress <= 1.0 ){

		std::cout << "[";
		int pos = barWidth * progress;
		for (int i = 0; i < barWidth; ++i) {
			if (i < pos) std::cout << "=";
			else if (i == pos) std::cout << ">";
			else std::cout << " ";
		}
		std::cout << "] " << int(progress * 100.0) << " %\r";
		std::cout.flush();
	}
}


std::vector< std::string > split( std::string s, char delim ){

	std::stringstream ssString( s );
	std::vector< std::string > splitString;
	std::string entry;

	while ( std::getline( ssString, entry, delim ) ){

		splitString.push_back( entry );
	}
	return splitString;
}


int argMin( std::vector< double > &vec ){

	double smallest = vec[0];
	int index = 0;

	for ( unsigned int i = 1; i < vec.size(); i++ ){

		if ( vec[i] < smallest ){

			smallest = vec[i];
			index = i;
		}
	}
	return index;
}
