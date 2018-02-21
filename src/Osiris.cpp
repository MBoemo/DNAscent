//----------------------------------------------------------
// Copyright 2017 University of Oxford
// Written by Michael A. Boemo (michael.boemo@path.ox.ac.uk)
// This software is licensed under GPL-3.0.  You should have
// received a copy of the license with this software.  If
// not, please Email the author.
//----------------------------------------------------------


#include <string>
#include <map>
#include <functional>
#include "Osiris_train.h"
#include "Osiris_fixedPos.h"
#include "Osiris_detect.h"
#include "error_handling.h"


/*prototype */
int show_options( int, char** );

/*map from name of the Osiris function passed as argument on the command line to the function that it should call */
static std::map< std::string, std::function< int( int, char** ) > > executables = {
	{"train", 	train_main},
	{"fixedPos", 	fixedPos_main},
	{"detect", 	detect_main},
	{"--help",	show_options},
	{"-h",	show_options}
};


int show_options( int, char** ){

	std::cout << "Osiris is a software tool for detecting base analogues in Oxford Nanopore reads.  To run Osiris, do: " \
	<< std::endl \
	<< "  ./Osiris [executable] [arguments]" \
	<< std::endl \
	<< "The executables that Osiris can run are:" \
	<< std::endl;

	for ( auto &exec : executables ){
		std::cout << "  "<< exec.first << std::endl;
	}

	return 0;

}


/*main Osiris executable that will link to other executables */
int main( int argc, char** argv ){

	if ( argc < 2 ){
		std::cout << "Exiting with error.  No Osiris executable specified." << std::endl <<  show_options( argc, argv );
		exit( EXIT_FAILURE );
	}

	std::string runThisExecutable( argv[ 1 ] );
	auto iter = executables.find( runThisExecutable );

	if ( iter == executables.end() ){
		std::cout << "Exiting with error.  Unknown Osiris executable specified." << std::endl <<  show_options( argc, argv );
		exit( EXIT_FAILURE );
	}

	return iter -> second( argc - 1, argv + 1 );

}
