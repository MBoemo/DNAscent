//----------------------------------------------------------
// Copyright 2017 University of Oxford
// Written by Michael A. Boemo (michael.boemo@path.ox.ac.uk)
//----------------------------------------------------------


#include "Osiris_fixedPos.h"
#include "common.h"
#include "build_model.h"
#include "data_IO.h"
#include "error_handling.h"


static const char *help=
"fixedPos: Osiris executable that trains for a base analogue's signal in a fixed position.\n"
"To run Osiris fixedPos, do:\n"
"  ./Osiris fixedPos [arguments]\n"
"Example:\n"
"  ./Osiris fixedPos -r /path/to/reference.fasta -p 29 -om /path/to/template_median68pA.model -d /path/to/data.foh -o output.txt -t 20 -sc 30\n"
"Required arguments are:\n"
"  -r,--reference            path to reference file in fasta format,\n"
"  -p,--position             position of analogue in training data,\n"
"  -om,--ont-model           path to 6mer pore model file (provided by ONT) over bases {A,T,G,C},\n"
"  -d,--trainingData         path to training data in the .foh format (can be made with Python Osiris),\n"
"  -o,--output               path to the training values at each position on the reference.\n"
"Optional arguments are:\n"
"  -t,--threads              number of threads (default is 1 thread),\n"
"  -sc,--soft-clipping       restrict training to this window size around a region of interest (default is off).\n";

struct Arguments {
	std::string referenceFilename;
	int analoguePosition;
	std::string trainingDataFilename;
	std::string ontModelFilename;
	std::string trainingOutputFilename;
	int threads;
	bool softClip;
	int SCwindow;
};

Arguments parseFixedPosArguments( int argc, char** argv ){

	if( argc < 2 ){

		std::cout << "Exiting with error.  Insufficient arguments passed to Osiris train." << std::endl << help << std::endl;
		exit(EXIT_FAILURE);

	}

	if ( std::string( argv[ 1 ] ) == "-h" or std::string( argv[ 1 ] ) == "--help" ){

		std::cout << help << std::endl;
		exit(EXIT_SUCCESS);
	}
	else if( argc < 4 ){

		std::cout << "Exiting with error.  Insufficient arguments passed to Osiris fixedPos." << std::endl;
		exit(EXIT_FAILURE);

	}

	Arguments trainArgs;

	/*defaults - we'll override these if the option was specified by the user */
	trainArgs.threads = 1;
	trainArgs.softClip = false;

	/*parse the command line arguments */
	for ( unsigned int i = 1; i < argc; ){

		std::string flag( argv[ i ] );

		if ( flag == "-r" or flag == "--reference" ){

			std::string strArg( argv[ i + 1 ] );
			trainArgs.referenceFilename = strArg;
			i+=2;	

		}
		else if ( flag == "-om" or flag == "--ont-model" ){

			std::string strArg( argv[ i + 1 ] );
			trainArgs.ontModelFilename = strArg;
			i+=2;

		}
		else if ( flag == "-d" or flag == "--trainingData" ){

			std::string strArg( argv[ i + 1 ] );
			trainArgs.trainingDataFilename = strArg;
			i+=2;

		}
		else if ( flag == "-o" or flag == "--output" ){

			std::string strArg( argv[ i + 1 ] );
			trainArgs.trainingOutputFilename = strArg;
			i+=2;

		}
		else if ( flag == "-t" or flag == "--threads" ){

			std::string strArg( argv[ i + 1 ] );
			trainArgs.threads = std::stoi( strArg.c_str() );
			i+=2;

		}
		else if ( flag == "-p" or flag == "--position" ){

			std::string strArg( argv[ i + 1 ] );
			trainArgs.analoguePosition = std::stoi( strArg.c_str() );
			i+=2;

		}
		else if ( flag == "-sc" or flag == "--soft-clipping" ){

			trainArgs.softClip = true;
			std::string strArg( argv[ i + 1 ] );
			trainArgs.SCwindow = std::stoi( strArg.c_str() );
			i+=2;

		}
		else throw InvalidOption( flag );
	}

	return trainArgs;

}


int fixedPos_main( int argc, char** argv ){

	Arguments trainArgs = parseFixedPosArguments( argc, argv );

	std::string reference = import_reference( trainArgs.referenceFilename );
	std::map< std::string, std::pair< double, double > > ontModel =  import_poreModel( trainArgs.ontModelFilename );
	std::map< std::string, std::vector< std::vector< double > > > trainingData = import_foh( trainArgs.trainingDataFilename );

	/*output file*/
	std::ofstream outFile;

	int brduDomLoc = trainArgs.analoguePosition;

	outFile.open( trainArgs.trainingOutputFilename );
	if ( not outFile.is_open() ) throw IOerror( trainArgs.trainingOutputFilename );

	for( auto iter = trainingData.cbegin(); iter != trainingData.cend(); ++iter ){

		std::vector< std::vector< double > > events = iter -> second;

		/*if soft clipping was specified, truncate the reference and events with dynamic time warping */
		if ( trainArgs.softClip == true ){

			if ( ( brduDomLoc - trainArgs.SCwindow < 0 ) or ( brduDomLoc + trainArgs.SCwindow > reference.length() ) ){

				std::cout << "Exiting with error.  Soft clipping window exceeds reference length.  Reduce window size." << std::endl;
				exit(EXIT_FAILURE);

			}

			reference = reference.substr( brduDomLoc - trainArgs.SCwindow, 6 + 2*trainArgs.SCwindow );
			brduDomLoc = trainArgs.SCwindow;
			events = filterEvents( reference, ontModel, events );
			
		}

		/*do the training */
		std::stringstream ss = buildAndTrainHMM( reference, ontModel, events, trainArgs.threads, true );

		/*if we specified that we want a log file, read the ss buffer into it now */
		std::stringstream ssLog( ss.str() );
		outFile << ssLog.rdbuf();

	}

	/*if we opened a log file to write on, close it now */
	outFile.close();

	return 0;

}
