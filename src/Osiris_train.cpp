//----------------------------------------------------------
// Copyright 2017 University of Oxford
// Written by Michael A. Boemo (michael.boemo@path.ox.ac.uk)
//----------------------------------------------------------


#include "Osiris_train.h"
#include "utility.h"


static const char *help=
"train: Osiris executable that determines the mean and standard deviation of a base analogue's current.\n"
"To run Osiris train, do:\n"
"  ./Osiris train [arguments]\n"
"Required arguments are:\n"
"  -r,--reference            full path to reference file in fasta format,\n"
"  -m,--model                full path to Oxford Nanopore pore model file,\n"
"  -d,--trainingData         full path to training data in the .foh format (can be made with Python Osiris\n"
"  -o,--output               full path to the output file which will contain the trained values.\n"
"Optional arguments are:\n"
"  -t,--threads              number of threads (default is 1 thread).\n";



struct Arguments {
	std::string referenceFilename;
	std::string trainingDataFilename;
	std::string baseModelFilename;
	std::string trainingOutputFilename;
	int threads;
};


Arguments parseArguments( int argc, char** argv ){

	if( argc < 2 ){

		std::cout << "Exiting with error.  Insufficient arguments passed to Osiris train." << std::endl << help << std::endl;
		exit(EXIT_FAILURE);

	}

	if ( std::string( argv[ 1 ] ) == "-h" or std::string( argv[ 1 ] ) == "--help" ){

		std::cout << help << std::endl;
		exit(EXIT_SUCCESS);
	}
	else if( argc < 4 ){

		std::cout << "Exiting with error.  Insufficient arguments passed to Osiris train." << std::endl;
		exit(EXIT_FAILURE);

	}

	Arguments trainArgs;

	/*defaults - we'll override these if the option was specified by the user */
	trainArgs.threads = 1;

	for ( unsigned int i = 1; i < argc; i+=2 ){

		std::string flag( argv[ i ] );

		if ( flag == "-r" or flag == "--reference" ){

			std::string strArg( argv[ i + 1 ] );
			trainArgs.referenceFilename = strArg;	

		}
		else if ( flag == "-m" or flag == "--model" ){

			std::string strArg( argv[ i + 1 ] );
			trainArgs.baseModelFilename = strArg;

		}
		else if ( flag == "-d" or flag == "--trainingData" ){

			std::string strArg( argv[ i + 1 ] );
			trainArgs.trainingDataFilename = strArg;

		}
		else if ( flag == "-o" or flag == "--output" ){

			std::string strArg( argv[ i + 1 ] );
			trainArgs.trainingOutputFilename = strArg;

		}
		else if ( flag == "-t" or flag == "--threads" ){

			std::string strArg( argv[ i + 1 ] );
			trainArgs.threads = std::stoi( strArg.c_str() );

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

	/*IO */
	std::ofstream outFile;
	outFile.open( trainArgs.trainingOutputFilename );
	if ( not outFile.is_open() ){
		std::cout << "Exiting with error.  Training data file could not be opened." << std::endl;
		exit(EXIT_FAILURE);
	}

	for( auto iter = trainingData.cbegin(); iter != trainingData.cend(); ++iter ){

		std::string refLocal = reference;

		std::string adenDomain = iter -> first;
		std::string brduDomain = reverseComplement( adenDomain );

		int adenDomLoc = refLocal.find( "NNNANNN" );
		int brduDomLoc = refLocal.find( "NNNTNNN" );

		refLocal.replace( adenDomLoc, adenDomain.length(), adenDomain );
		refLocal.replace( brduDomLoc, brduDomain.length(), brduDomain );

		std::vector< std::vector< double > > events = iter -> second;

		std::stringstream ss = buildAndTrainHMM( refLocal, baseModel, events, trainArgs.threads );

		outFile << ">" << adenDomain << std::endl << ss.rdbuf();

	}

	outFile.close();

	return 0;

}
