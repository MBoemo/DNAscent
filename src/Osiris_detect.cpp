//----------------------------------------------------------
// Copyright 2017 University of Oxford
// Written by Michael A. Boemo (michael.boemo@path.ox.ac.uk)
//----------------------------------------------------------


#include "Osiris_detect.h"
#include "common.h"
#include "build_model.h"
#include "data_IO.h"


static const char *help=
"detect: Osiris executable that determines the position of base analogues in a Nanopore read.\n"
"To run Osiris detect, do:\n"
"  ./Osiris detect [arguments]\n"
"Example:\n"
"  ./Osiris detect -bm /path/to/template_median68pA.model -am /path/to/analogue.model -d /path/to/data.foh -o output.txt -t 20\n"
"Required arguments are:\n"
"  -om,--ont-model           full path to 6mer pore model file (provided by ONT) over bases {A,T,G,C},\n"
"  -am,--analogue-model      full path to 6mer pore model file that includes analogues,\n"
"  -d,--data                 full path to .fdh file to run detection on,\n"
"  -o,--output               full path to the output file which will contain the calls.\n"
"Optional arguments are:\n"
"  -t,--threads              number of threads (default is 1 thread).\n";

struct Arguments {
	std::string ontModelFilename;
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
	for ( unsigned int i = 1; i < argc; ){

		std::string flag( argv[ i ] );

		if ( flag == "-om" or flag == "--ont-model" ){

			std::string strArg( argv[ i + 1 ] );
			trainArgs.ontModelFilename = strArg;
			i+=2;

		}
		else if ( flag == "-am" or flag == "--analogue-model" ){

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
		else{

			std::cout << "Exiting with error.  Invalid argument passed to Osiris train." << std::endl;
			exit(EXIT_FAILURE);

		}

	}

	return trainArgs;

}


int detect_main( int argc, char** argv ){

	Arguments trainArgs = parseDetectionArguments( argc, argv );

	std::map< std::string, std::pair< double, double > > ontModel =  import_poreModel( trainArgs.ontModelFilename );
	std::map< std::string, std::pair< double, double > > analogueModel =  import_poreModel( trainArgs.analogueModelFilename );
	std::vector< detectionTuple > detectData = import_fdh( trainArgs.dataFilename );

	/*IO */
	std::ofstream outFile;
	outFile.open( trainArgs.outputFilename );
	if ( not outFile.is_open() ){
		std::cout << "Exiting with error.  Output file could not be opened." << std::endl;
		exit(EXIT_FAILURE);
	}

	std::stringstream ss; //out stringstream

	for( unsigned int i = 0; i < detectData.size(); i++ ){

		/*take the read, and if it's really short, toss it out */
		std::string basecalls = detectData[ i ].basecalls;
		if ( basecalls.length() < 40 ){
			continue;
		}
		std::string fast5File = detectData[ i ].filename;
		std::vector< double > events = detectData[ i ].events;

		ss << ">" << fast5File << std::endl;

		/*make an ideal signal from the basecalls */
		std::vector< double > generatedSignal = generateSignal( basecalls, ontModel );

		/*get a rough alignment from the basecalls to events with dynamic time warping */
		std::map< int, std::vector< int > > warpPath = dynamicTimewarping( events, generatedSignal );

		/*if we can't get an alignment that goes to the beginning of the read, don't try to detect anything in this read and skip to the next one */
		if ( warpPath.count( 10 ) == 0 ){
			continue;
		}

		/*bound these - things tend to be rubbish at the start and end so start 15 bases in and end 30 bases early */
		for ( unsigned int j = 20; j < basecalls.length() - 30; j++ ){
			
			if ( basecalls.substr( j, 1 ) == "T"){

				std::string readSnippet = basecalls.substr( j - 10, 26 );
				std::vector< double >::const_iterator first = events.begin() + warpPath[ j - 10 ].back();
				std::vector< double >::const_iterator last = events.begin() + warpPath[ j + 16 ].front();				
				std::vector< double > eventSnippet( first, last );

				double logProbThymidine = buildAndDetectHMM( readSnippet, ontModel, analogueModel, eventSnippet, false );
				double logProbAnalogue = buildAndDetectHMM( readSnippet, ontModel, analogueModel, eventSnippet, true );

				double logLikelihoodRatio = logProbAnalogue - logProbThymidine;

				ss << j << "\t" << logLikelihoodRatio << std::endl;

			}

		}

		/*write the log probabilities for this file to the output file */
		outFile << ss.rdbuf();

	}

	outFile.close();

	return 0;

}
