//----------------------------------------------------------
// Copyright 2017 University of Oxford
// Written by Michael A. Boemo (michael.boemo@path.ox.ac.uk)
// This software is licensed under GPL-3.0.  You should have
// received a copy of the license with this software.  If
// not, please Email the author.
//----------------------------------------------------------

#include <exception>
#include <math.h>
#include <iostream>
#include <fstream>
#include "common.h"
#include "data_IO.h"
#include "error_handling.h"
#include "event_handling.h"
#include "../Penthus/src/error_handling.h"
#include "../Penthus/src/hmm.h"
#include "../Penthus/src/states.h"
#include "poreModels.h"
#include "Osiris_scan.h"
#include "poreSpecificParameters.h"

static const char *help=
"train: Osiris executable that scans a sequence for events that deviate from expectation.\n"
"To run Osiris scan, do:\n"
"  ./Osiris scan [arguments]\n"
"Example:\n"
"  ./Osiris scan -d /path/to/data.foh -o output.txt -t 20\n"
"Required arguments are:\n"
"  -d,--trainingData         path to training data in the .foh format (made with prepTrainingData.py),\n"
"  -o,--output               path to the output pore model file that Osiris will train.\n"
"Optional arguments are:\n"
"  -t,--threads              number of threads (default is 1 thread).\n";

struct Arguments {

	std::string trainingDataFilename;
	std::string trainingOutputFilename;
	bool logFile;
	std::string logFilename;
	int threads;
};

Arguments parseTrainingArguments( int argc, char** argv ){

	if( argc < 2 ) throw InsufficientArguments();

	if ( std::string( argv[ 1 ] ) == "-h" or std::string( argv[ 1 ] ) == "--help" ){

		std::cout << help << std::endl;
		exit(EXIT_SUCCESS);
	}
	//else if( argc < 4 ) throw InsufficientArguments();

	Arguments trainArgs;

	/*defaults - we'll override these if the option was specified by the user */
	trainArgs.threads = 1;

	/*parse the command line arguments */
	for ( int i = 1; i < argc; ){

		std::string flag( argv[ i ] );

		if ( flag == "-d" or flag == "--trainingData" ){

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
		else throw InvalidOption( flag );
	}
	return trainArgs;
}


int train_main( int argc, char** argv ){

	Arguments trainArgs = parseTrainingArguments( argc, argv );

	/*get a filestream to the foh file - we'll load training data dynamically */
	std::ifstream fohFile( trainArgs.trainingDataFilename );
	if ( not fohFile.is_open() ) throw IOerror( trainArgs.trainingDataFilename );

	/*open output file */
	std::ofstream outFile( trainArgs.trainingOutputFilename );
	if ( not fohFile.is_open() ) throw IOerror( trainArgs.trainingOutputFilename );

	/*read the foh header - total count and reference */
	std::string line, reference;
	std::getline( fohFile, reference );
	std::getline( fohFile, line );
	int trainingTotal = atoi(line.c_str());


	int prog = 0;

	while ( std::getline( fohFile, line) ){

		/*get data for a read from foh */
		read currentRead;

		/*the basecall line */
		currentRead.basecalls = line;

		/*the bounds line */
		std::getline( fohFile, line );
		(currentRead.ROIbounds).first = atoi( (line.substr( 0, line.find(' ') )).c_str() );
		(currentRead.ROIbounds).second = atoi( (line.substr( line.find(' ') + 1, line.size() - line.find(' ') )).c_str() );

		/*the raw signal line */
		std::getline( fohFile, line );
		std::vector< double > rawSignals;
		std::istringstream ss( line );
		std::string event;
		while ( std::getline( ss, event, ' ' ) ){

			rawSignals.push_back( atof( event.c_str() ) );
		}
		currentRead.raw = rawSignals;

		/*normalise for shift and scale */
		eventDataForRead eventData = normaliseEvents( currentRead );

		/*disregard this event if the quality score is too low */
		if ( fabs(eventData.qualityScore) > 1.0 ) continue;

		std::string refSeqMapped = reference.substr((currentRead.ROIbounds).first, (currentRead.ROIbounds).second - (currentRead.ROIbounds).first);

		HiddenMarkovModel hmm = HiddenMarkovModel( 3*refSeqMapped.length(), 3*refSeqMapped.length() + 2 );

		std::pair< double, double > emissionMeanAndStd;

		/*STATES - vector (of vectors) to hold the states at each position on the reference - fill with dummy values */
		std::vector< std::vector< State > > states( 6, std::vector< State >( refSeqMapped.length() - 5, State( NULL, "", "", "", 1.0 ) ) );

		/*DISTRIBUTIONS - vector to hold normal distributions, a single uniform and silent distribution to use for everything else */
		std::vector< NormalDistribution > nd;
		nd.reserve( refSeqMapped.length() - 5 );		
		SilentDistribution sd( 0.0, 0.0 );
		UniformDistribution ud( 50.0, 150.0 );

		std::string loc;

		/*create the distributions that we need */			
		for ( unsigned int i = 0; i < refSeqMapped.length() - 5; i++ ){

			emissionMeanAndStd = FiveMer_model.at( refSeqMapped.substr( i, 5 ) );		
			nd.push_back( NormalDistribution( emissionMeanAndStd.first, emissionMeanAndStd.second ) );		
		}

		/*add states to model, handle internal module transitions */
		for ( unsigned int i = 0; i < refSeqMapped.length() - 5; i++ ){

			loc = std::to_string( i + (currentRead.ROIbounds).first );
			std::string fiveMer = refSeqMapped.substr( i, 5 );

			states[ 0 ][ i ] = State( &sd, 		loc + "_SS",	fiveMer,	"", 		1.0 );
			states[ 1 ][ i ] = State( &sd,		loc + "_D", 	fiveMer,	"", 		1.0 );		
			states[ 2 ][ i ] = State( &ud,		loc + "_I", 	fiveMer,	"", 		1.0 );
			states[ 3 ][ i ] = State( &nd[i], 	loc + "_M1", 	fiveMer,	loc + "_match", 1.0 );
			states[ 4 ][ i ] = State( &nd[i], 	loc + "_M2", 	fiveMer,	loc + "_match", 1.0 );
			states[ 5 ][ i ] = State( &sd, 		loc + "_SE", 	fiveMer,	"", 		1.0 );

			/*add state to the model */
			for ( unsigned int j = 0; j < 6; j++ ){

				states[ j ][ i ].meta = fiveMer;
				hmm.add_state( states[ j ][ i ] );
			}

			/*transitions between states, internal to a single base */
			/*from SS */
			hmm.add_transition( states[0][i], states[3][i], internalSS2M1 );
			hmm.add_transition( states[0][i], states[4][i], internalSS2M2 );

			/*from D */
			hmm.add_transition( states[1][i], states[2][i], internalD2I );

			/*from I */
			hmm.add_transition( states[2][i], states[2][i], internalI2I );
			hmm.add_transition( states[2][i], states[0][i], internalI2SS );

			/*from M1 */
			hmm.add_transition( states[3][i], states[3][i], internalM12M1 );
			hmm.add_transition( states[3][i], states[5][i], internalM12SE );

			/*from M2 */
			hmm.add_transition( states[4][i], states[4][i], internalM22M2 );
			hmm.add_transition( states[4][i], states[5][i], internalM22SE );

			/*from SE */
			hmm.add_transition( states[5][i], states[2][i], internalSE2I );		

		}

		/*add transitions between modules (external transitions) */
		for ( unsigned int i = 0; i < refSeqMapped.length() - 6; i++ ){

			hmm.add_transition( states[1][i], states[1][i + 1], externalD2D );
			hmm.add_transition( states[1][i], states[0][i + 1], externalD2SS );
			hmm.add_transition( states[2][i], states[0][i + 1], externalI2SS );
			hmm.add_transition( states[5][i], states[0][i + 1], externalSE2SS );
			hmm.add_transition( states[5][i], states[1][i + 1], externalSE2D );
		}

		/*handle start states */
		hmm.add_transition( hmm.start, states[0][0], 0.5 );
		hmm.add_transition( hmm.start, states[1][0], 0.5 );

		/*handle end states */
		hmm.add_transition( states[1][refSeqMapped.length() - 6], hmm.end, externalD2D + externalD2SS );
		hmm.add_transition( states[2][refSeqMapped.length() - 6], hmm.end, externalI2SS );
		hmm.add_transition( states[5][refSeqMapped.length() - 6], hmm.end, externalSE2SS + externalSE2D );

		hmm.finalise();

		/*event alignment using Viterbi */
		auto viterbiData = hmm.viterbi( eventData.normalisedEvents );
		std::vector< State > statePath = viterbiData.second;

		/*filter state path to only emitting states */
		std::vector< State > emittingStatePath;
		for ( auto s = statePath.begin(); s < statePath.end(); s++ ){

			/*check character */
			if ( ((*s).name).substr(((*s).name).find('_') + 1,1) == "M" or ((*s).name).substr(((*s).name).find('_') + 1,1) == "I" ){
				emittingStatePath.push_back( *s );
			}
		}

		/*print out the alignment to output */
		for ( unsigned int i = 0; i < (eventData.normalisedEvents).size(); i++ ){

			outFile << (eventData.normalisedEvents)[i] << '\t' << emittingStatePath[i].name << '\t' << emittingStatePath[i].meta << '\t' << emittingStatePath[i].dist -> param1 << '\t' << emittingStatePath[i].dist -> param2 << std::endl;
		}
		outFile << ">" << std::endl;

		displayProgress( prog, trainingTotal );
		prog++;
	}
	return 0;
}
