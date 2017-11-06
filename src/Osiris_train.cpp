//----------------------------------------------------------
// Copyright 2017 University of Oxford
// Written by Michael A. Boemo (michael.boemo@path.ox.ac.uk)
// This software is licensed under GPL-3.0.  You should have
// received a copy of the license with this software.  If
// not, please Email the author.
//----------------------------------------------------------


#include <exception>
#include "Osiris_train.h"
#include "common.h"
#include "build_model.h"
#include "data_IO.h"
#include "error_handling.h"
#include "event_handling.h"
#include "../Penthus/src/error_handling.h"

static const char *help=
"train: Osiris executable that determines the mean and standard deviation of a base analogue's current.\n"
"To run Osiris train, do:\n"
"  ./Osiris train [arguments]\n"
"Example:\n"
"  ./Osiris train -r /path/to/reference.fasta -p 3and4 -bm /path/to/template_median68pA.model -d /path/to/data.foh -o output.txt -t 20 -sc 30\n"
"Required arguments are:\n"
"  -r,--reference            path to reference file in fasta format,\n"
"  -p,--position             position of analogue in training data (valid arguments are 1and2, 3and4, or 5and6),\n"
"  -d,--trainingData         path to training data in the .foh format (made with prepHairpinData.py),\n"
"  -o,--output               path to the output pore model file that Osiris will train.\n"
"Optional arguments are:\n"
"  -t,--threads              number of threads (default is 1 thread),\n"
"  -sc,--soft-clipping       restrict training to this window size around a region of interest (default is off),\n"
"  -l,--log-file             training log file for the training values at each position on the reference (default is none).\n";

struct Arguments {
	std::string referenceFilename;
	std::string analoguePosition;
	std::string trainingDataFilename;
	std::string trainingOutputFilename;
	bool logFile;
	std::string logFilename;
	int threads;
	bool softClip;
	int SCwindow;
};

Arguments parseTrainingArguments( int argc, char** argv ){

	if( argc < 2 ) throw InsufficientArguments();

	if ( std::string( argv[ 1 ] ) == "-h" or std::string( argv[ 1 ] ) == "--help" ){

		std::cout << help << std::endl;
		exit(EXIT_SUCCESS);
	}
	else if( argc < 4 ) throw InsufficientArguments();

	Arguments trainArgs;

	/*defaults - we'll override these if the option was specified by the user */
	trainArgs.threads = 1;
	trainArgs.softClip = false;
	trainArgs.logFile = false;

	/*parse the command line arguments */
	for ( int i = 1; i < argc; ){

		std::string flag( argv[ i ] );

		if ( flag == "-r" or flag == "--reference" ){

			std::string strArg( argv[ i + 1 ] );
			trainArgs.referenceFilename = strArg;
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
			trainArgs.analoguePosition = strArg;
			i+=2;
		}
		else if ( flag == "-sc" or flag == "--soft-clipping" ){

			trainArgs.softClip = true;
			std::string strArg( argv[ i + 1 ] );
			trainArgs.SCwindow = std::stoi( strArg.c_str() );
			i+=2;
		}
		else if ( flag == "-l" or flag == "--log-file" ){

			trainArgs.logFile = true;
			std::string strArg( argv[ i + 1 ] );
			trainArgs.logFilename = strArg;
			i+=2;
		}
		else throw InvalidOption( flag );
	}
	return trainArgs;
}


int train_main( int argc, char** argv ){

	Arguments trainArgs = parseTrainingArguments( argc, argv );

	std::string reference = import_reference( trainArgs.referenceFilename );
	//std::map< std::string, std::vector< std::vector< double > > > trainingData = import_foh( trainArgs.trainingDataFilename );

	/*MODIFIED HERE */
	std::map< std::string, std::vector< std::vector< double > > > normalisedEvents = segmentEvents( trainArgs.trainingDataFilename );

	/*END */


	return 0;

}
