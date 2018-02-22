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
"  -m,--analogue-model       path to 5mer pore model file that includes analogues,\n"
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

	/*import the analogue pore model that we specified on the command line */
	std::map< std::string, std::pair< double, double > > analogueModel =  import_poreModel( trainArgs.analogueModelFilename );

	/*IO */
	std::ofstream outFile;
	outFile.open( trainArgs.outputFilename );
	if ( not outFile.is_open() ) throw IOerror( trainArgs.outputFilename );
	std::ifstream fdhFile( trainArgs.dataFilename );
	if ( not fdhFile.is_open() ) throw IOerror( trainArgs.dataFilename );
	std::stringstream ss;

	/*read foh header */
	std::string line;
	std::getline( fdhFile, line );
	std::string strDetectionTotal = line.substr(0,line.find(' '));
	int detectionTotal = atoi(strDetectionTotal.c_str());

	unsigned int windowLength = 20;
	int prog = 0;
	std::vector< read > buffer;
	std::string detectGroup;
	std::cout << "Detecting base analogues..." << std::endl;

	/*get a single detection read from the input fdh */
	while ( std::getline( fdhFile, detectGroup, '<' ) ){

		assert( buffer.size() < trainArgs.threads );

		/*push that read to the buffer */
		if ( not fdhFile.eof() ) buffer.push_back( getDetectionFrom_fdh( detectGroup ) );

		/*if we've filled up the buffer, run detection on the number of threads we specified */
		if ( buffer.size() == trainArgs.threads or fdhFile.eof() ){

			#pragma omp parallel for default(none) schedule(dynamic) shared(detectionTotal, prog, windowLength, FiveMer_model, analogueModel, trainArgs, buffer, outFile) private(ss) num_threads(trainArgs.threads)
			for( auto r = buffer.begin(); r < buffer.end(); r++ ){

				/*normalise raw for shift and scale, and get an alignment between 5mers and the events that produced them */
				eventDataForRead eventData = normaliseEvents( *r, false );
				
				/*disregard this read if we don't meet the minimum required quality score for the 5mer alignment */				
				if ( eventData.qualityScore > 4.0 ) continue;

				/*push the filename for this read to the output */
				ss << ">" << r -> filename << std::endl;

				int readHead = 0;

				/*exclude the starts and ends of the read, as the alignment tends to be worse there */
				for ( unsigned int i = 2*windowLength; i < (r -> basecalls).length() - 2*windowLength; i++ ){
			
					/*for each T we find in the read, calculate the log likelihood that it's an analogue */
					if ( (r -> basecalls).substr( i, 1 ) == "T"){

						bool makeCall = true;

						/*make sure we have all the 5mers we need */
						std::vector< unsigned int > analoguePositions = {i, i-1, i-2, i-3};
						for ( auto pos = analoguePositions.begin(); pos < analoguePositions.end(); pos++){

							if (analogueModel.count(((r -> basecalls).substr(*pos, 5)).replace(i - *pos, 1, "B")) == 0) makeCall = false;
						}

						/*if we do have all the 5mers that we need, then calculate the log likelihood of an analogue at this position */
						if ( makeCall ){
							std::string readSnippet = (r -> basecalls).substr(i - windowLength, 2*windowLength);
							std::vector< double > eventSnippet;

							/*get the events that correspond to the read snippet */
							for ( unsigned int j = readHead; j < (eventData.eventAlignment).size(); j++ ){

								/*move the readhead so that we can navigate through events efficiently */
								if ( (eventData.eventAlignment)[j].second == i - windowLength ) readHead = j;

								/*if an event has been aligned to a position in the window, add it */
								if ( (eventData.eventAlignment)[j].second >= i - windowLength and (eventData.eventAlignment)[j].second <= i + windowLength ){

									eventSnippet.push_back( (eventData.normalisedEvents)[j] );
								}

								/*stop once we get to the end of the window */
								if ( (eventData.eventAlignment)[j].second > i + windowLength ) break;
							}
							double logProbThymidine = buildAndDetectHMM( readSnippet, FiveMer_model, analogueModel, eventSnippet, windowLength, false );
							double logProbAnalogue = buildAndDetectHMM( readSnippet, FiveMer_model, analogueModel, eventSnippet, windowLength, true );
							double logLikelihoodRatio = logProbAnalogue - logProbThymidine;

							ss << i << "\t" << logLikelihoodRatio << std::endl;
						}
					}
				}

				#pragma omp atomic 
				prog++;

				#pragma omp critical
				{	/*write the log probabilities for this file to the output file */
					outFile << ss.rdbuf();
					displayProgress( prog, detectionTotal );
				}
			}

			/*empty the buffer */
			buffer.clear();
		}
		if ( fdhFile.eof() ) break;

	}
	displayProgress( detectionTotal, detectionTotal );
	std::cout << "\nDone." << std::endl;

	outFile.close();
	fdhFile.close();

	return 0;
}
