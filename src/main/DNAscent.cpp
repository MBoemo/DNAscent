//----------------------------------------------------------
// Copyright 2019 University of Oxford
// Written by Michael A. Boemo (mb915@cam.ac.uk)
// This software is licensed under GPL-3.0.  You should have
// received a copy of the license with this software.  If
// not, please Email the author.
//----------------------------------------------------------

#include <string>
#include <iostream>
#include <map>
#include <functional>
#include "../detect.h"
#include "../forkSense.h"
#include "../index.h"
#include "../common.h"
#include "../trainCNN.h"
#include "../alignment.h"
#include "../trainGMM.h"
#include "../config.h"

/*prototype */
int show_options_DNAscent( int, char** );

/*map from name of the DNAscent function passed as argument on the command line to the function that it should call */
static std::map< std::string, std::function< int( int, char** ) > > executables = {
	{"index", 	index_main},
	{"detect", 	detect_main},
	{"forkSense", 	sense_main},
	{"align", 	align_main},
	{"trainCNN", 	data_main},
	{"trainGMM", 	train_main},
	{"--help",	show_options_DNAscent},
	{"-h",		show_options_DNAscent},
	{"-v",		show_version},
	{"--version",	show_version}
};


int show_options_DNAscent( int, char** ){

	std::cout << "DNAscent is a software tool for detecting regions of base analogue incorporation in Oxford Nanopore reads.\nTo run DNAscent, do: " \
	<< std::endl \
	<< "  ./DNAscent [executable] [arguments]" \
	<< std::endl \
	<< "The executables that DNAscent can run are:" \
	<< std::endl;

	for ( auto &exec : executables ){
		std::cout << "  "<< exec.first << std::endl;
	}
	std::cout << "Version: " << VERSION << std::endl;
	std::cout << "Written by Michael Boemo, Department of Pathology, University of Cambridge." << std::endl;
	std::cout << "Please submit bug reports to GitHub Issues." << std::endl;
	return 0;
}

Global_Config Pore_Substrate_Config;

/*main DNAscent executable that will link to other executables */
int main( int argc, char** argv ){

	//configure for DNA R10
	Pore_Substrate_Config.configure_DNA_R10();

	if ( argc < 2 ){
		std::cout << "Exiting with error.  No DNAscent executable specified." << std::endl <<  show_options_DNAscent( argc, argv );
		exit( EXIT_FAILURE );
	}

	std::string runThisExecutable( argv[ 1 ] );
	auto iter = executables.find( runThisExecutable );

	if ( iter == executables.end() ){
		std::cout << "Exiting with error.  Unknown DNAscent executable specified." << std::endl <<  show_options_DNAscent( argc, argv );
		exit( EXIT_FAILURE );
	}

	return iter -> second( argc - 1, argv + 1 );
}