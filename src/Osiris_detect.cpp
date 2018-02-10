//----------------------------------------------------------
// Copyright 2017 University of Oxford
// Written by Michael A. Boemo (michael.boemo@path.ox.ac.uk)
// This software is licensed under GPL-3.0.  You should have
// received a copy of the license with this software.  If
// not, please Email the author.
//----------------------------------------------------------


#include <fstream>
#include "Osiris_detect.h"
#include "common.h"
#include "build_model.h"
#include "data_IO.h"
#include "error_handling.h"
#include "event_handling.h"
#include "poreModels.h"


static const char *help=
"detect: Osiris executable that determines the position of base analogues in a Nanopore read.\n"
"To run Osiris detect, do:\n"
"  ./Osiris detect [arguments]\n"
"Example:\n"
"  ./Osiris detect -m /path/to/analogue.model -d /path/to/data.foh -o output.txt -t 20\n"
"Required arguments are:\n"
"  -m,--analogue-model       path to 6mer pore model file that includes analogues,\n"
"  -d,--data                 path to .fdh file to run detection on,\n"
"  -o,--output               path to the output file which will contain the calls.\n"
"Optional arguments are:\n"
"  -t,--threads              number of threads (default is 1 thread).\n";

struct Arguments {
	std::string analogueModelFilename;
	std::string dataFilename;
	std::string outputFilename;
	int threads;
};

Arguments parseDetectionArguments( int argc, char** argv ){

	if( argc < 2 ){

		std::cout << "Exiting with error.  Insufficient arguments passed to Osiris train." << std::endl << help << std::endl;
		exit(EXIT_FAILURE);

	}

	if ( std::string( argv[ 1 ] ) == "-h" or std::string( argv[ 1 ] ) == "--help" ){

		std::cout << help << std::endl;
		exit(EXIT_SUCCESS);
	}
	else if( argc < 4 ){

		std::cout << "Exiting with error.  Insufficient arguments passed to Osiris detect." << std::endl;
		exit(EXIT_FAILURE);

	}

	Arguments trainArgs;

	/*defaults - we'll override these if the option was specified by the user */
	trainArgs.threads = 1;

	/*parse the command line arguments */
	for ( int i = 1; i < argc; ){

		std::string flag( argv[ i ] );

		if ( flag == "-m" or flag == "--analogue-model" ){

			std::string strArg( argv[ i + 1 ] );
			trainArgs.analogueModelFilename = strArg;
			i+=2;

		}
		else if ( flag == "-o" or flag == "--output" ){

			std::string strArg( argv[ i + 1 ] );
			trainArgs.outputFilename = strArg;
			i+=2;

		}
		else if ( flag == "-t" or flag == "--threads" ){

			std::string strArg( argv[ i + 1 ] );
			trainArgs.threads = std::stoi( strArg.c_str() );
			i+=2;

		}
		else if ( flag == "-d" or flag == "--data" ){

			std::string strArg( argv[ i + 1 ] );
			trainArgs.dataFilename = strArg;
			i+=2;

		}
		else throw InvalidOption( flag );
	}

	return trainArgs;

}


int detect_main( int argc, char** argv ){

	Arguments trainArgs = parseDetectionArguments( argc, argv );

	std::map< std::string, std::pair< double, double > > analogueModel =  import_poreModel( trainArgs.analogueModelFilename );
	std::vector< read > detectData = import_fdh( trainArgs.dataFilename );

	/*IO */
	std::ofstream outFile;
	outFile.open( trainArgs.outputFilename );
	if ( not outFile.is_open() ) throw IOerror( trainArgs.outputFilename );

	std::stringstream ss; //out stringstream

	unsigned int windowLength = 10;
	int prog = 0;
	int total = detectData.size();
	std::cout << "Detecting base analogues..." << std::endl;

	#pragma omp parallel for default(none) schedule(dynamic) shared(total, prog, windowLength, FiveMer_model, analogueModel, trainArgs, detectData, outFile) private(ss) num_threads(trainArgs.threads)
	for( auto r = detectData.begin(); r < detectData.end(); r++ ){

		/*normalise raw for shift and scale, and get an alignment between 5mers and the events that produced them */
		eventDataForRead eventData = normaliseEvents( *r, false );

		int readHead = 0;

		/*bound these - things tend to be rubbish at the start and end so start 15 bases in and end 30 bases early */
		for ( unsigned int i = windowLength; i < (r -> basecalls).length() - 3*windowLength; i++ ){
			
			if ( (r -> basecalls).substr( i, 1 ) == "T"){

				std::string readSnippet = (r -> basecalls).substr( i - windowLength, 2*windowLength + 5 );
				std::vector< double > eventSnippet;

				/*get the events that correspond to the read snippet */
				for ( unsigned int j = readHead; j < (eventData.eventAlignment).size(); j++ ){

					if ( (eventData.eventAlignment)[j].second == i - windowLength ){

						readHead = j;
					}

					if ( (eventData.eventAlignment)[j].second >= i - windowLength and (eventData.eventAlignment)[j].second <= i + windowLength + 5 ){

						eventSnippet.push_back( (eventData.normalisedEvents)[j] );
					}

					if ( (eventData.eventAlignment)[j].second > i + windowLength + 5 ) break;
				}


				double logProbThymidine = buildAndDetectHMM( readSnippet, FiveMer_model, analogueModel, eventSnippet, false );
				double logProbAnalogue = buildAndDetectHMM( readSnippet, FiveMer_model, analogueModel, eventSnippet, true );

				double logLikelihoodRatio = logProbAnalogue - logProbThymidine;

				ss << i << "\t" << logLikelihoodRatio << std::endl;

			}
		}

		#pragma omp atomic 
		prog++;

		#pragma omp critical
		{	/*write the log probabilities for this file to the output file */
			outFile << ss.rdbuf();
			displayProgress( prog, total );
		}
	}
	displayProgress( prog, total );
	std::cout << "\nDone." << std::endl;

	outFile.close();

	return 0;

}
