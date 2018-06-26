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
#include "../Penthus/src/unsupervised_learning.h"
#include "poreModels.h"
#include "Osiris_train.h"
#include "poreSpecificParameters.h"

static const char *help=
"train: Osiris executable that determines the mean and standard deviation of a base analogue's current.\n"
"To run Osiris train, do:\n"
"  ./Osiris train [arguments]\n"
"Example:\n"
"  ./Osiris train -d /path/to/data.foh -o output.txt -t 20\n"
"Required arguments are:\n"
"  -d,--trainingData         path to training data in the .foh format (made with prepTrainingData.py),\n"
"  -o,--output               path to the output pore model file that Osiris will train.\n"
"Optional arguments are:\n"
"  -e,--eventalign           take input from nanopolish eventalign,\n"
"  -t,--threads              number of threads (default is 1 thread).\n";

struct Arguments {

	std::string trainingDataFilename;
	std::string trainingOutputFilename;
	bool logFile;
	std::string logFilename;
	int threads;
	bool useNanopolish;
	std::string eventalignFilename;
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
	trainArgs.useNanopolish = false;

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
		else if ( flag == "-e" or flag == "--eventalign" ){

			std::string strArg( argv[ i + 1 ] );
			trainArgs.eventalignFilename = strArg;
			trainArgs.useNanopolish = true;
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


void alignToSixmer( Arguments &trainArgs ){

	/*get a filestream to the foh file - we'll load training data dynamically */
	std::ifstream fohFile( trainArgs.trainingDataFilename );
	if ( not fohFile.is_open() ) throw IOerror( trainArgs.trainingDataFilename );

	/*read the foh header - total count */
	std::string line;
	std::getline( fohFile, line );
	int trainingTotal = atoi(line.c_str());

	/*initialise progress */
	progressBar pb_align( trainingTotal );
	int prog = 0, offloadCount = 0, failed = 0;

	/*buffers */
	std::map< std::string, std::vector< double > > eventPileup;
	std::vector< read > buffer;

	/*open work file */
	std::ofstream workFile( "workingData.osiris" );
	if ( not workFile.is_open() ) throw IOerror( "workingData.osiris" );

	/*align the training data */
	std::cout << "Aligning events..." << std::endl;
	while ( std::getline( fohFile, line) ){
		
		/*get data for a read from foh */
		read currentRead;

		/*subsequence of the reference that the read mapped to */
		currentRead.mappedRefSubseq = line;

		/*the basecalled sequence */
		std::getline( fohFile, line );
		currentRead.basecall = line;

		/*the raw signal line */
		std::getline( fohFile, line );
		std::vector< double > rawSignals;
		std::istringstream ss( line );
		std::string event;
		while ( std::getline( ss, event, ' ' ) ){

			rawSignals.push_back( atof( event.c_str() ) );
		}
		currentRead.raw = rawSignals;

		/*push it to the buffer or run Viterbi if the buffer is full */
		buffer.push_back( currentRead );
		if ( (buffer.size() < trainArgs.threads)  ) continue;

		/*Viterbi event alignment */
		#pragma omp parallel for schedule(dynamic) shared(failed,pb_align,eventPileup,buffer,trainArgs,trainingTotal,prog,SixMer_model,internalM12I, internalI2I, internalM12M1, externalD2D, externalD2M1, externalI2M1, externalM12D, externalM12M1) num_threads(trainArgs.threads)
		for ( auto r = buffer.begin(); r < buffer.end(); r++ ){

			/*normalise for shift and scale */
			eventDataForRead eventData = normaliseEvents( *r );

			/*disregard this event if the quality score is too low */
			//if ( fabs(eventData.qualityScore) > 1.0 ){

			//	failed++;
			//	prog++;
			//	continue;
			//}

			/*get the subsequence of the reference this read mapped to and build an HMM from it */
			std::string refSeqMapped = (*r).mappedRefSubseq;
			HiddenMarkovModel hmm = HiddenMarkovModel( 3*refSeqMapped.length(), 3*refSeqMapped.length() + 2 );

			/*STATES - vector (of vectors) to hold the states at each position on the reference - fill with dummy values */
			std::vector< std::vector< State > > states( 6, std::vector< State >( refSeqMapped.length() - 5, State( NULL, "", "", "", 1.0 ) ) );

			/*DISTRIBUTIONS - vector to hold normal distributions, a single uniform and silent distribution to use for everything else */
			std::vector< NormalDistribution > nd;
			nd.reserve( refSeqMapped.length() - 5 );

			SilentDistribution sd( 0.0, 0.0 );
			UniformDistribution ud( 50.0, 150.0 );

			std::string loc, sixMer;

			/*create make normal distributions for each reference position using the ONT 6mer model */			
			for ( unsigned int i = 0; i < refSeqMapped.length() - 6; i++ ){

				sixMer = refSeqMapped.substr( i, 6 );
				nd.push_back( NormalDistribution( SixMer_model[sixMer].first, SixMer_model[sixMer].second ) );
			}

			/*the first insertion state after start */
			State firstI = State( &ud, "-1_I", "", "", 1.0 );
			hmm.add_state( firstI );

			/*add states to the model, handle internal module transitions */
			for ( unsigned int i = 0; i < refSeqMapped.length() - 6; i++ ){

				loc = std::to_string( i );
				sixMer = refSeqMapped.substr( i, 6 );

				states[ 0 ][ i ] = State( &sd,		loc + "_D", 	sixMer,	"", 		1.0 );		
				states[ 1 ][ i ] = State( &ud,		loc + "_I", 	sixMer,	"", 		1.0 );
				states[ 2 ][ i ] = State( &nd[i], 	loc + "_M1", 	sixMer,	loc + "_match", 1.0 );

				/*add state to the model */
				for ( unsigned int j = 0; j < 3; j++ ){

					states[ j ][ i ].meta = sixMer;
					hmm.add_state( states[ j ][ i ] );
				}

				/*transitions between states, internal to a single base */
				/*from I */
				hmm.add_transition( states[1][i], states[1][i], internalI2I );

				/*from M1 */
				hmm.add_transition( states[2][i], states[2][i], internalM12M1 );
				hmm.add_transition( states[2][i], states[1][i], internalM12I );
			}

			/*add transitions between modules (external transitions) */
			for ( unsigned int i = 0; i < refSeqMapped.length() - 7; i++ ){

				/*from D */
				hmm.add_transition( states[0][i], states[0][i + 1], externalD2D );
				hmm.add_transition( states[0][i], states[2][i + 1], externalD2M1 );

				/*from I */
				hmm.add_transition( states[1][i], states[2][i + 1], externalI2M1 );

				/*from M */
				hmm.add_transition( states[2][i], states[0][i + 1], externalM12D );
				hmm.add_transition( states[2][i], states[2][i + 1], externalM12M1 );
			}

			/*handle start states */
			hmm.add_transition( hmm.start, firstI, 0.25 );
			hmm.add_transition( hmm.start, states[0][0], 0.25 );
			hmm.add_transition( hmm.start, states[2][0], 0.5 );

			/*transitions from first insertion */
			hmm.add_transition( firstI, firstI, 0.25 );
			hmm.add_transition( firstI, states[0][0], 0.25 );
			hmm.add_transition( firstI, states[2][0], 0.5 );

			/*handle end states */
			hmm.add_transition( states[0][refSeqMapped.length() - 7], hmm.end, 1.0 );
			hmm.add_transition( states[1][refSeqMapped.length() - 7], hmm.end, externalI2M1 );
			hmm.add_transition( states[2][refSeqMapped.length() - 7], hmm.end, externalM12M1 + externalM12D );

			hmm.finalise();

			/*do the event alignment with the Penthus Viterbi algorithm */
			auto viterbiData = hmm.viterbi( eventData.normalisedEvents ); 
			double viterbiScore = viterbiData.first;
			if  ( std::isnan( viterbiScore ) ){

				failed++;
				prog++;
				continue;
			}	
			std::vector< std::string > statePath = viterbiData.second;

			/*filter state path to only emitting states */
			std::vector< std::string > emittingStatePath;
			for ( auto s = statePath.begin(); s < statePath.end(); s++ ){

				if ( (*s).substr((*s).find('_') + 1,1) == "M" or (*s).substr((*s).find('_') + 1,1) == "I" ){
					emittingStatePath.push_back( *s );
				}
			}

			/*add the events to the eventPileup hash table - these events are keyed by their position in the reference */
			#pragma omp critical
			{
				for ( unsigned int i = 0; i < (eventData.normalisedEvents).size(); i++ ){

					std::vector< std::string > findIndex = split( emittingStatePath[i], '_' );
					unsigned int posOnReference = atoi(findIndex[0].c_str());

					if ( findIndex[1].substr(0,1) == "M" ){
						//std::cout << (eventData.normalisedEvents)[i] << '\t' << posOnReference << '\t' << emittingStatePath[i] << '\t' << reference.substr(posOnReference,5) << '\t' << FiveMer_model[reference.substr(posOnReference,5)].first << '\t' << FiveMer_model[reference.substr(posOnReference,5)].second << std::endl;
						eventPileup[refSeqMapped.substr(posOnReference,6)].push_back( (eventData.normalisedEvents)[i] );
					}
				}
			pb_align.displayProgress( prog, failed );
			prog++;
			}
		}
		offloadCount++;
		buffer.clear();

		/*after running through the buffer a few times, offload event data to file so we don't flood memory */
		if ( offloadCount == 5 ){

			for ( auto iter = eventPileup.cbegin(); iter != eventPileup.cend(); ++iter ){

				workFile << iter -> first;
				for ( auto e = (iter -> second).begin(); e < (iter -> second).end(); e++ ){

					workFile << ' ' << *e;
				}
				workFile << std::endl;
			}
			offloadCount = 0;
			eventPileup.clear();
		}
	}

	/*empty the buffer at the end */
	for ( auto iter = eventPileup.cbegin(); iter != eventPileup.cend(); ++iter ){

		workFile << iter -> first;
		for ( auto e = (iter -> second).begin(); e < (iter -> second).end(); e++ ){

			workFile << ' ' << *e;
		}
		workFile << std::endl;
	}
	eventPileup.clear();

	/*wrap up */
	fohFile.close();
	workFile.close();
	eventPileup.clear();
	pb_align.displayProgress( trainingTotal, failed );
	failed = prog = 0;
	std::cout << std::endl << "Done." << std::endl;
}


int train_main( int argc, char** argv ){

	Arguments trainArgs = parseTrainingArguments( argc, argv );

	if ( not trainArgs.useNanopolish){

		alignToSixmer( trainArgs );
	}

	/*open output file */
	std::ofstream outFile( trainArgs.trainingOutputFilename );
	if ( not outFile.is_open() ) throw IOerror( trainArgs.trainingOutputFilename );

	/*get events from the work file */
	std::cout << "Fitting Gaussian mixture model..." << std::endl;

	std::string line;
	int prog, failed;

	/*fudge for openmp */
	std::map< int, std::string > indexToSixmer;
	std::map< std::string, int > sixmerToIndex;
	int index = 0;
	for ( auto i = SixMer_model.cbegin(); i != SixMer_model.cend(); i++ ){

		indexToSixmer[index] = i -> first;
		sixmerToIndex[i->first] = index;
		index++;
	}

	std::vector< std::vector< double > > importedEvents( 4096 );

	if (trainArgs.useNanopolish) {

		std::ifstream eventFile(trainArgs.eventalignFilename);
		if ( not eventFile.is_open() ) throw IOerror( trainArgs.eventalignFilename );

		std::getline( eventFile, line);//throw away the header

		while ( std::getline( eventFile, line) ){

			std::istringstream ss( line );
			std::string sixMer, entry;
			double eventMean, eventLength;

			int col = 0;
			while ( std::getline( ss, entry, '\t' ) ){

				if ( col == 9 ){

					sixMer = entry;
					break;
				}
				else if ( col == 6 ){

					eventMean = atof( entry.c_str() );
				}
				else if ( col == 8 ){

					eventLength = atof( entry.c_str() );
				}
				col++;
			}
			if ( eventLength >= 0.002 ) importedEvents[sixmerToIndex[sixMer]].push_back( eventMean );
		}
		eventFile.close();
	}
	else {

		std::ifstream eventFile("workingData.osiris");
		if ( not eventFile.is_open() ) throw IOerror( "workingData.osiris" );

		while ( std::getline( eventFile, line) ){

			std::vector< double > rawSignals;
			std::istringstream ss( line );
			std::string event;
			std::getline( ss, event, ' ' );
			std::string sixMer = event;

			while ( std::getline( ss, event, ' ' ) ){

				rawSignals.push_back( atof( event.c_str() ) );
			}
			importedEvents[sixmerToIndex[sixMer]].insert( importedEvents[sixmerToIndex[sixMer]].end(), rawSignals.begin(), rawSignals.end() );
		}
		eventFile.close();
	}

	/*fit a mixture model to the events that aligned to each position in the reference */
	outFile << "6mer" << '\t' << "ONT_mean" << '\t' << "ONT_stdv" << '\t' << "pi_1" << '\t' << "mean_1" << '\t' << "stdv_1" << '\t' << "pi_2" << '\t' << "mean_2" << '\t' << "stdv_2" << std::endl;
	progressBar pb_fit( importedEvents.size() );

	#pragma omp parallel for schedule(dynamic) shared(pb_fit, indexToSixmer, SixMer_model, prog, failed, outFile, importedEvents, trainArgs) num_threads(trainArgs.threads)
	for ( int i = 0; i < importedEvents.size(); i++ ){

		/*don't train if we have less than 2000 events for this 6mer */
		if ( importedEvents[i].size() < 200 ){

			prog++;
			continue;
		}

		std::string sixMer = indexToSixmer[i];
		double mu1, stdv1, mu2, stdv2;//, mu3, stdv3;

		/*get the ONT distribution for the mixture */
		mu1 = SixMer_model[sixMer].first;
		stdv1 = SixMer_model[sixMer].second;

		/*make a second distribution that's similar to the ONT distribution */
		mu2 = SixMer_model[sixMer].first;
		stdv2 = 2*SixMer_model[sixMer].second;

		/*make a second distribution that's similar to the ONT distribution */
		//mu3 = SixMer_model[sixMer].first;
		//stdv3 = 5*SixMer_model[sixMer].second;

		/*fit the model */
		std::vector< double > fitParameters;
		try{
			//fitParameters = gaussianMixtureEM_TRIMODAL( mu1, stdv1, mu2, stdv2, mu3, stdv3, importedEvents[i], 0.0001 );
			fitParameters = gaussianMixtureEM_PRIOR( mu1, stdv1, mu2, stdv2, importedEvents[i], 10 );
		}
		catch ( NegativeLog &nl ){

			failed++;
			prog++;
			continue;
		}
		#pragma omp critical
		{	
			//outFile << sixMer << '\t' << SixMer_model[sixMer].first << '\t' << SixMer_model[sixMer].second << '\t' << fitParameters[0] << '\t' << fitParameters[1] << '\t' << fitParameters[2] << '\t' << fitParameters[3] << '\t' << fitParameters[4] << '\t' << fitParameters[5] << '\t' << fitParameters[6] << '\t' << fitParameters[7] << '\t' << fitParameters[8] << '\t' << (importedEvents[i]).size() << std::endl;
			outFile << sixMer << '\t' << SixMer_model[sixMer].first << '\t' << SixMer_model[sixMer].second << '\t' << fitParameters[0] << '\t' << fitParameters[1] << '\t' << fitParameters[2] << '\t' << fitParameters[3] << '\t' << fitParameters[4] << '\t' << fitParameters[5] << '\t' << (importedEvents[i]).size() << std::endl; 

			pb_fit.displayProgress( prog, failed );
		}
		prog++;
	}
	outFile.close();
	std::cout << std::endl << "Done." << std::endl;

	return 0;
}
