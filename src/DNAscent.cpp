//----------------------------------------------------------
// Copyright 2019 University of Oxford
// Written by Michael A. Boemo (michael.boemo@path.ox.ac.uk)
// This software is licensed under GPL-2.0.  You should have
// received a copy of the license with this software.  If
// not, please Email the author.
//----------------------------------------------------------

#include <string>
#include <iostream>
#include <map>
#include <functional>
#include "detect.h"
#include "regions.h"
#include "psl.h"
#include "index.h"


/*prototype */
int show_options( int, char** );

/*map from name of the DNAscent function passed as argument on the command line to the function that it should call */
static std::map< std::string, std::function< int( int, char** ) > > executables = {
	{"detect", 	detect_main},
	{"psl", 	psl_main},
	{"regions", 	regions_main},
	{"index", 	index_main},
	{"--help",	show_options},
	{"-h",	show_options}
};


int show_options( int, char** ){

	std::cout << "DNAscent is a software tool for detecting regions of base analogue incorporation in Oxford Nanopore reads.\nTo run DNAscent, do: " \
	<< std::endl \
	<< "  ./DNAscent [executable] [arguments]" \
	<< std::endl \
	<< "The executables that DNAscent can run are:" \
	<< std::endl;

	for ( auto &exec : executables ){
		std::cout << "  "<< exec.first << std::endl;
	}

	return 0;

}


/*main DNAscent executable that will link to other executables */
int main( int argc, char** argv ){

	if ( argc < 2 ){
		std::cout << "Exiting with error.  No DNAscent executable specified." << std::endl <<  show_options( argc, argv );
		exit( EXIT_FAILURE );
	}

	std::string runThisExecutable( argv[ 1 ] );
	auto iter = executables.find( runThisExecutable );

	if ( iter == executables.end() ){
		std::cout << "Exiting with error.  Unknown DNAscent executable specified." << std::endl <<  show_options( argc, argv );
		exit( EXIT_FAILURE );
	}

	return iter -> second( argc - 1, argv + 1 );
}
