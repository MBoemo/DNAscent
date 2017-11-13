//----------------------------------------------------------
// Copyright 2017 University of Oxford
// Written by Michael A. Boemo (michael.boemo@path.ox.ac.uk)
// This software is licensed under GPL-3.0.  You should have
// received a copy of the license with this software.  If
// not, please Email the author.
//----------------------------------------------------------

#include <exception>
#include <fstream>
#include "common.h"
#include "build_model.h"
#include "data_IO.h"
#include "error_handling.h"
#include "event_handling.h"
#include "../Penthus/src/error_handling.h"
#include "poreModels.h"
#include "Osiris_train.h"


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
"  -c,--clipping             restrict training to a small window around the region of interest,\n"
"  -l,--log-file             training log file for the training values at each position on the reference (default is none).\n";

struct Arguments {
	std::string referenceFilename;
	std::string analoguePosition;
	std::string trainingDataFilename;
	std::string trainingOutputFilename;
	bool logFile;
	std::string logFilename;
	int threads;
	bool clip;
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
	trainArgs.clip = false;
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
		else if ( flag == "-c" or flag == "--clipping" ){

			trainArgs.clip = true;
			i++;
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

	/*import a reference from a fasta file and training data from a foh file.  Normalise the training data. */
	std::string reference = import_reference( trainArgs.referenceFilename );
	std::map< std::string, std::vector< std::vector< double > > > trainingData = segmentEvents( trainArgs.trainingDataFilename, trainArgs.threads, trainArgs.clip );

	/*log file IO - if we specified that we want a log file, open it here */
	std::ofstream logFile;
	if ( trainArgs.logFile == true ){
		logFile.open( trainArgs.logFilename );
		if ( not logFile.is_open() ) throw IOerror( trainArgs.logFilename );
	}

	std::map< std::string, std::vector< double> > trainedModel; //this will be filled up with results from Penthus
	int prog = 0;

	/*iterate on the 7mers that we want to train on */
	for( auto iter = trainingData.cbegin(); iter != trainingData.cend(); ++iter ){

		std::string refLocal = reference;
		std::vector< std::vector< double > > events = iter -> second;

		std::string brduDomain = iter -> first;
		std::string adenDomain = reverseComplement( brduDomain );

		/*our reference has N's in it: replace them with the appropriate 7mer that we're goin to train on */
		int positionNorm, adenDomLoc, brduDomLoc;

		if ( trainArgs.analoguePosition == "1and2" ){
			adenDomLoc = refLocal.find( "NNNNNAN" );
			brduDomLoc = refLocal.find( "NTNNNNN" );
			positionNorm = 0;
		}
		else if ( trainArgs.analoguePosition == "3and4" ){
			adenDomLoc = refLocal.find( "NNNANNN" );
			brduDomLoc = refLocal.find( "NNNTNNN" );
			positionNorm = 2;
		}
		else if ( trainArgs.analoguePosition == "5and6" ){
			adenDomLoc = refLocal.find( "NANNNNN" );
			brduDomLoc = refLocal.find( "NNNNNTN" );
			positionNorm = 4;
		}
		else{
			std::cout << "Exiting with error.  Invalid option passed with the -p or --position flag.  Valid options are 1and2, 3and4, or 5and6." << std::endl;
			exit(EXIT_FAILURE);
		}
	
		refLocal.replace( adenDomLoc, adenDomain.length(), adenDomain );
		refLocal.replace( brduDomLoc, brduDomain.length(), brduDomain );
		
		/*check for unresolved Ns and fail if there are any */
		for ( auto it = refLocal.begin(); it < refLocal.end(); it++ ){
			if ( *it == 'N' ){
				std::cout << "Exiting with error.  Reference contains unresolved N's." << std::endl;
				exit(EXIT_FAILURE);
			}
		}

		/*if clipping was specified, adjust the reference accordingly */
		if ( trainArgs.clip == true ){

			assert( ( brduDomLoc - 6 >= 0 ) and ( brduDomLoc + 13 <= refLocal.length() ) );

			refLocal = refLocal.substr( brduDomLoc - 6, 19 );
			brduDomLoc = 6;
		}

		/*do the training */
		std::stringstream ss;
		try{
			ss = buildAndTrainHMM( refLocal, SixMer_model, events, trainArgs.threads, false );
		}
		catch ( NumericalInstability &ni ){
			std::cout << ni.what() << std::endl << "Aborted training on this 6mer, skipping: "<< brduDomain << std::endl;
			displayProgress( prog, trainingData.size() );
			prog++;
			continue; 
		}

		/*if we specified that we want a log file, read the ss buffer into it now */
		std::stringstream ssLog( ss.str() );
		if ( trainArgs.logFile == true ){
			logFile << ">" << brduDomain << std::endl << ssLog.rdbuf();
		}

		/*hacky bodge to get the training data out at the relevant position without making Penthus specialised */
		std::string line;
		while ( std::getline( ss, line ) ){
			
			std::vector< std::string > findIndex = split( line, '_' );
			unsigned int i = atoi(findIndex[0].c_str());

			if ( i == brduDomLoc ){

				std::string atFirstPos = ( brduDomain.substr( 0, 6 ) ).replace( positionNorm + 1, 1, "B" );
				std::vector< std::string > splitLine = split( line, '\t' );
				
				if ( trainedModel.count( atFirstPos ) > 0 ){

					if ( atof( splitLine[ 5 ].c_str() ) < trainedModel[ atFirstPos ][ 1 ] ){

						trainedModel[ atFirstPos ] = { atof( splitLine[ 3 ].c_str() ), atof( splitLine[ 5 ].c_str() ), atof( splitLine[ 2 ].c_str() ), atof( splitLine[ 4 ].c_str() ) };
					}
				}
				else{

					trainedModel[ atFirstPos ] = { atof( splitLine[ 3 ].c_str() ), atof( splitLine[ 5 ].c_str() ), atof( splitLine[ 2 ].c_str() ), atof( splitLine[ 4 ].c_str() ) };
				}

			}
			else if ( i == brduDomLoc + 1 ){

				std::string atSecondPos = ( brduDomain.substr( 1, 6 ) ).replace( positionNorm, 1, "B" );
				std::vector< std::string > splitLine = split( line, '\t' );

				if ( trainedModel.count( atSecondPos ) > 0 ){

					if ( atof( splitLine[ 5 ].c_str() ) < trainedModel[ atSecondPos ][ 1 ] ){

						trainedModel[ atSecondPos ] = { atof( splitLine[ 3 ].c_str() ), atof( splitLine[ 5 ].c_str() ), atof( splitLine[ 2 ].c_str() ), atof( splitLine[ 4 ].c_str() ) };
					}
				}
				else{

					trainedModel[ atSecondPos ] = { atof( splitLine[ 3 ].c_str() ), atof( splitLine[ 5 ].c_str() ), atof( splitLine[ 2 ].c_str() ), atof( splitLine[ 4 ].c_str() ) };
				}
			}
			else if ( i > brduDomLoc + 1 ){
		
				break;
			}
		}
		displayProgress( prog, trainingData.size() );
		prog++;
	}

	/*make a pore model file from the map */
	export_poreModel( trainedModel, trainArgs.trainingOutputFilename );

	/*if we opened a log file to write on, close it now */
	if ( trainArgs.logFile == true ){
		logFile.close();
	}

	/*some wrap-up messages */
	std::cout << std::endl << "Done." << std::endl;


	return 0;

}
