//----------------------------------------------------------
// Copyright 2017 University of Oxford
// Written by Michael A. Boemo (michael.boemo@path.ox.ac.uk)
//----------------------------------------------------------


#include "Osiris_train.h"


struct Arguments {
	std::string ref;
	//training data
	std::string baseModelFilename;
};


Arguments parseArguments( int argc, char** argv ){

	if ( argc < 4 ){

		std::cout << "Exiting with error.  Insufficient arguments passed to Osiris train." << std::endl;
		exit(EXIT_FAILURE);

	}

	Arguments trainArgs;

	for ( unsigned int i = 1; i < argc; i+=2 ){

		std::string flag( argv[ i ] );

		if ( flag == "-r" ){

			std::string strArg( argv[ i + 1 ] );
			trainArgs.ref = strArg;	

		}
		else if ( flag == "-m" ){

			std::string strArg( argv[ i + 1 ] );
			trainArgs.baseModelFilename = strArg;

		}
		else{

			std::cout << "Exiting with error.  Invalid argument passed to Osiris train." << std::endl;
			exit(EXIT_FAILURE);

		}

	}

	return trainArgs;

}


int train_main( int argc, char** argv ){

	Arguments trainArgs = parseArguments( argc, argv );

	std::string reference = import_reference( trainArgs.ref );
	std::map< std::string, std::pair< double, double > > baseModel =  import_poreModel( trainArgs.baseModelFilename );

	HiddenMarkovModel hmm = build_trainingHMM( reference, baseModel );

	return 0;

}
