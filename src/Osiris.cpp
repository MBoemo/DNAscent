//----------------------------------------------------------
// Copyright 2017 University of Oxford
// Written by Michael A. Boemo (michael.boemo@path.ox.ac.uk)
//----------------------------------------------------------


#include <string>
#include <map>
#include <functional>

#include "Osiris_train.h"

/*map from name of the Osiris function passed as argument on the command line to the function that it should call */
static std::map< std::string, std::function< int( int, char** ) > > executables = {
	{"train", 	train_main}
};


/*main Osiris executable that will link to other executables */
int main( int argc, char** argv ){

	std::string runThisExecutable( argv[ 1 ] );
	auto iter = executables.find( runThisExecutable );

	return iter -> second( argc - 1, argv + 1 );

}
