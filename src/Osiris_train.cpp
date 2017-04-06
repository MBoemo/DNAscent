//----------------------------------------------------------
// Copyright 2017 University of Oxford
// Written by Michael A. Boemo (michael.boemo@path.ox.ac.uk)
//----------------------------------------------------------


#include "Osiris_train.h"


struct Arguments {
	std::string referenceFilename;
	std::string trainingDataFilename;
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
			trainArgs.referenceFilename = strArg;	

		}
		else if ( flag == "-m" ){

			std::string strArg( argv[ i + 1 ] );
			trainArgs.baseModelFilename = strArg;

		}
		else if ( flag == "-d" ){

			std::string strArg( argv[ i + 1 ] );
			trainArgs.trainingDataFilename = strArg;

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

	std::string reference = import_reference( trainArgs.referenceFilename );
	std::map< std::string, std::pair< double, double > > baseModel =  import_poreModel( trainArgs.baseModelFilename );

	std::map< std::string, std::vector< std::vector< double > > > trainingData = import_foh( trainArgs.trainingDataFilename );

	HiddenMarkovModel hmm;

	for( auto iter = trainingData.cbegin(); iter != trainingData.cend(); ++iter ){

		/*write something that replaces the reference NNNTNNN and NNNANNN with appropriate values from the key */

		hmm = build_trainingHMM( reference, baseModel );
		
		hmm.BaumWelch( iter -> second, 1.0 );

	}

	return 0;

}
