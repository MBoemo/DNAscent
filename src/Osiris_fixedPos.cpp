//----------------------------------------------------------
// Copyright 2017 University of Oxford
// Written by Michael A. Boemo (michael.boemo@path.ox.ac.uk)
// This software is licensed under GPL-3.0.  You should have
// received a copy of the license with this software.  If
// not, please Email the author.
//----------------------------------------------------------


#include <fstream>
#include "Osiris_fixedPos.h"
#include "common.h"
#include "build_model.h"
#include "data_IO.h"
#include "poreModels.h"
#include "error_handling.h"
#include "event_handling.h"
#include "../Penthus/src/error_handling.h"


static const char *help=
"fixedPos: Osiris executable that trains for a base analogue's signal in a fixed position.\n"
"To run Osiris fixedPos, do:\n"
"  ./Osiris fixedPos [arguments]\n"
"Example:\n"
"  ./Osiris fixedPos -r /path/to/reference.fasta -p 29 -om /path/to/template_median68pA.model -d /path/to/data.foh -o output.txt -t 20 -sc 30\n"
"Required arguments are:\n"
"  -r,--reference            path to reference file in fasta format,\n"
"  -p,--position             position of analogue in training data,\n"
"  -d,--trainingData         path to training data in the .foh format (can be made with Python Osiris),\n"
"  -o,--output               path to the training values at each position on the reference.\n"
"Optional arguments are:\n"
"  -t,--threads              number of threads (default is 1 thread),\n"
"  -c,--clipping             restrict training to a small window around the region of interest.\n"
"  -l,--log-file             training log file for the training values at each position on the reference (default is none).\n";

struct Arguments {
	std::string referenceFilename;
	unsigned int analoguePosition;
	std::string trainingDataFilename;
	std::string trainingOutputFilename;
	bool logFile;
	std::string logFilename;
	int threads;
	bool clip;
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
	trainArgs.clip = false;

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
			trainArgs.analoguePosition = std::stoi( strArg.c_str() );
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


int fixedPos_main( int argc, char** argv ){

	Arguments trainArgs = parseFixedPosArguments( argc, argv );

	/*import a reference from the input fasta file */
	std::string reference = import_reference( trainArgs.referenceFilename );

	/*get a filestream to the foh file - because this is fixed position, we don't care about the foh header */
	std::ifstream fohFile( trainArgs.trainingDataFilename );
	if ( not fohFile.is_open() ) throw IOerror( trainArgs.trainingDataFilename );

	/*log file IO - if we specified that we want a log file, open it here */
	std::ofstream logFile;
	if ( trainArgs.logFile == true ){
		logFile.open( trainArgs.logFilename );
		if ( not logFile.is_open() ) throw IOerror( trainArgs.logFilename );
	}

	/*we're looking for BrdU at a fixed position, at the position specified on the command line */
	unsigned int brduDomLoc = trainArgs.analoguePosition;

	/*this will be filled up with results from Penthus */
	std::map< std::string, std::vector< double> > trainedModel;

	/*iterate on the 7mers that we want to train on */
	std::pair< std::string, std::vector< read > > trainingGroup;
	std::string trainingGroupStr;
	while ( std::getline( fohFile, trainingGroupStr, '<' ) ){

		/*stop if we're at the end of the file */
		if ( fohFile.eof() ) break; 

		/*grab the training data from the foh */
		trainingGroup = getTrainingFrom_foh( trainingGroupStr );

		/*for each read, segment the events, normalise for shift and scale, and clip if specified */
		std::vector< std::vector< double > > events;
		bool clip = trainArgs.clip;
		int threads = trainArgs.threads;
		#pragma omp parallel for default(none) shared(trainingGroup, events, clip) num_threads(threads)
		for ( auto r = (trainingGroup.second).begin(); r < (trainingGroup.second).end(); r++ ){
			
			eventDataForRead localEvents = normaliseEvents( *r, clip );			
			
			#pragma omp critical
			events.push_back( localEvents.normalisedEvents );
		}

		/*if clipping was specified, adjust the reference accordingly */
		if ( trainArgs.clip == true ){

			/*make sure the domain is far enough from the 5' end that we have enough room to clip */
			assert( ( brduDomLoc - 16 >= 0 ) and ( brduDomLoc + 23 <= reference.length() ) );

			/*cut down the reference and reset brduDomLoc appropriately */
			reference = reference.substr( brduDomLoc - 16, 39 );
			brduDomLoc = 16;
		}

		/*do the training - we just have one reference so use the HMM verbose output */
		std::stringstream ss;
		try{
			ss = buildAndTrainHMM( reference, FiveMer_model, events, trainArgs.threads, true );
		}
		catch ( NumericalInstability &ni ){
			std::cout << ni.what() << std::endl << "Aborted training due to numerical instability - try adding more reads." << std::endl;
			continue; 
		}

		/*if we specified that we want a log file, read the ss buffer into it now */
		std::stringstream ssLog( ss.str() );
		if ( trainArgs.logFile == true ){
			logFile << ">fixedPositionTrained" << std::endl << ssLog.rdbuf();
		}

		/*use brduDomLoc to determine the positions we should look at to write the model */
		std::vector< unsigned int > posToLookAt = {brduDomLoc, brduDomLoc - 1, brduDomLoc - 2, brduDomLoc - 3, brduDomLoc - 4};

		/*bodge to get the training data out at the relevant position without making Penthus specialised */
		std::string line;
		while ( std::getline( ss, line ) ){
			
			/*get the position */
			std::vector< std::string > findIndex = split( line, '_' );
			unsigned int i = atoi(findIndex[0].c_str());

			if ( std::find(posToLookAt.begin(), posToLookAt.end(), i) != posToLookAt.end() ){

				std::vector< std::string > splitLine = split( line, '\t' );
				std::string fiveMer = reference.substr(i, 5);
				trainedModel[fiveMer] = {atof(splitLine[3].c_str()), atof(splitLine[5].c_str()), atof(splitLine[2].c_str()), atof(splitLine[4].c_str())};
			}
			else if ( i > posToLookAt.back() ) break;
		}
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
