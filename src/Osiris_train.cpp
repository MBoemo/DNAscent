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
"  ./Osiris train -r /path/to/reference.fasta -p 234 -d /path/to/data.foh -o output.txt -t 20 -c -l /path/to/myLogFile.log\n"
"Required arguments are:\n"
"  -r,--reference            path to reference file in fasta format,\n"
"  -p,--position             position of analogue in training data (valid arguments are 12, 234, or 45),\n"
"  -d,--trainingData         path to training data in the .foh format (made with prepTrainingData.py),\n"
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

	/*import a reference from the fasta file specified */
	std::string reference = import_reference( trainArgs.referenceFilename );

	/*get a filestream to the foh file - we'll load training data dynamically */
	std::ifstream fohFile( trainArgs.trainingDataFilename );
	if ( not fohFile.is_open() ) throw IOerror( trainArgs.trainingDataFilename );

	/*read the foh header */
	std::string line;
	std::getline( fohFile, line );
	std::string strTrainingTotal = line.substr(0,line.find(' '));
	int trainingTotal = atoi(strTrainingTotal.c_str());

	/*log file IO - if we specified that we want a log file, open it here */
	std::ofstream logFile;
	if ( trainArgs.logFile == true ){
		logFile.open( trainArgs.logFilename );
		if ( not logFile.is_open() ) throw IOerror( trainArgs.logFilename );
	}

	std::map< std::string, std::vector< double> > trainedModel; //this will be filled up with results from Penthus
	int prog = 0;

	std::cout << "Training..." << std::endl;

	/*iterate on the 7mers that we want to train on */
	std::pair< std::string, std::vector< read > > trainingGroup;
	std::string trainingGroupStr;
	while ( std::getline( fohFile, trainingGroupStr, '<' ) ){

		/*stop if we're at the end of the file */
		if ( fohFile.eof() ) break; 

		/*grab a 7mer's worth of training data from the foh */
		trainingGroup = getTrainingFrom_foh( trainingGroupStr );
		std::string brduDomain = trainingGroup.first;

		/*for each read, segment the events, normalise for shift and scale, and clip if specified */
		std::vector< std::vector< double > > events;
		bool clip = trainArgs.clip;
		int threads = trainArgs.threads;
		#pragma omp parallel for default(none) shared(trainingGroup, events, clip) num_threads(threads)
		for ( auto r = (trainingGroup.second).begin(); r < (trainingGroup.second).end(); r++ ){
			
			std::vector< double > localEvents = normaliseEvents( *r, clip );			
			
			#pragma omp critical
			events.push_back( localEvents );
		}

		/*our reference has N's in it: replace them with the appropriate 7mer that we're goin to train on */
		std::string refLocal = reference;
		unsigned int brduDomLoc;
		std::string brduDom_B_replace_T = brduDomain;

		std::vector< unsigned int > posToLookAt;
		if ( trainArgs.analoguePosition == "12" ){

			brduDomLoc = refLocal.find( "NTNNNNN" );
			posToLookAt = {0, 1};
			brduDom_B_replace_T.replace(1,1,"B");
		}
		else if ( trainArgs.analoguePosition == "234" ){

			brduDomLoc = refLocal.find( "NNNTNNN" );
			posToLookAt = {0, 1, 2};
			brduDom_B_replace_T.replace(3,1,"B");
		}
		else if ( trainArgs.analoguePosition == "45" ){

			brduDomLoc = refLocal.find( "NNNNNTN" );
			posToLookAt = {1, 2};
			brduDom_B_replace_T.replace(5,1,"B");
		}
		else{
			std::cout << "Exiting with error.  Invalid option passed with the -p or --position flag.  Valid options are 1and2, 3and4, or 5and6." << std::endl;
			exit(EXIT_FAILURE);
		}
	
		refLocal.replace( brduDomLoc, brduDomain.length(), brduDomain );
		
		/*check for unresolved Ns and fail if there are any */
		if ( refLocal.find("N") != std::string::npos ) throw NsInReference();

		/*if clipping was specified, adjust the reference accordingly */
		if ( trainArgs.clip == true ){

			/*make sure the domain is far enough from the 5' end that we have enough room to clip */
			assert( ( brduDomLoc - 16 >= 0 ) and ( brduDomLoc + 23 <= refLocal.length() ) );

			/*cut down the reference and reset brduDomLoc appropriately */
			refLocal = refLocal.substr( brduDomLoc - 16, 39 );
			brduDomLoc = 16;
		}

		/*fix posToLookAt with brduDomLoc that we set */
		for ( unsigned int i = 0; i < posToLookAt.size(); i++ ){

			posToLookAt[i] = posToLookAt[i] + brduDomLoc;
		}

		/*do the training */
		std::stringstream ss;
		try{
			ss = buildAndTrainHMM( refLocal, FiveMer_model, events, trainArgs.threads, false );
		}
		catch ( NumericalInstability &ni ){
			std::cout << ni.what() << std::endl << "Aborted training on this 7mer, skipping: "<< brduDom_B_replace_T << std::endl;
			displayProgress( prog, trainingTotal );
			prog++;
			continue; 
		}

		/*if we specified that we want a log file, read the ss buffer into it now */
		std::stringstream ssLog( ss.str() );
		if ( trainArgs.logFile == true ){
			logFile << ">" << brduDom_B_replace_T << std::endl << ssLog.rdbuf();
		}

		/*bodge to get the training data out at the relevant position without making Penthus specialised */
		std::string line;
		while ( std::getline( ss, line ) ){
			
			/*get the position */
			std::vector< std::string > findIndex = split( line, '_' );
			unsigned int i = atoi(findIndex[0].c_str());

			if ( std::find(posToLookAt.begin(), posToLookAt.end(), i) != posToLookAt.end() ){

				std::vector< std::string > splitLine = split( line, '\t' );
				std::string fiveMer = brduDom_B_replace_T.substr(i-brduDomLoc, 5);
				if ( trainedModel.count( fiveMer ) > 0 ){

					if ( atof(splitLine[ 5 ].c_str()) < trainedModel[fiveMer][1] ){

						trainedModel[fiveMer] = {atof(splitLine[3].c_str()), atof(splitLine[5].c_str()), atof(splitLine[2].c_str()), atof(splitLine[4].c_str())};
					}
				}
				else{

					trainedModel[fiveMer] = {atof(splitLine[3].c_str()), atof(splitLine[5].c_str()), atof(splitLine[2].c_str()), atof(splitLine[4].c_str())};
				}
			}
			else if ( i > posToLookAt.back() ) break;
		}
		displayProgress( prog, trainingTotal );
		prog++;
	}
	displayProgress( trainingTotal, trainingTotal );

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
