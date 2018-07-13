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
#include "poreSpecificParameters.h"


static const char *help=
"detect: Osiris executable that determines the position of base analogues in a Nanopore read.\n"
"To run Osiris detect, do:\n"
"  ./Osiris detect [arguments]\n"
"Example:\n"
"  ./Osiris detect -m /path/to/analogue.model -d /path/to/data.foh -o output.txt -t 20\n"
"Required arguments are:\n"
"  -m,--analogue-model       path to 5mer pore model file that includes analogues,\n"
"  -d,--data                 path to .foh file to run detection on,\n"
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


double seqProbability( std::string &sequence, std::vector<double> &events, std::map< std::string, std::pair< double, double > > &analogueModel, int BrdU_idx, bool isBrdU ){


	HiddenMarkovModel hmm = HiddenMarkovModel( 3*sequence.length(), 3*sequence.length() + 2 );

	/*STATES - vector (of vectors) to hold the states at each position on the reference - fill with dummy values */
	std::vector< std::vector< State > > states( 6, std::vector< State >( sequence.length() - 5, State( NULL, "", "", "", 1.0 ) ) );

	/*DISTRIBUTIONS - vector to hold normal distributions, a single uniform and silent distribution to use for everything else */
	std::vector< NormalDistribution > nd;
	nd.reserve( sequence.length() - 5 );

	SilentDistribution sd( 0.0, 0.0 );
	UniformDistribution ud( 50.0, 150.0 );

	std::string loc, sixMer;
		
	/*create make normal distributions for each reference position using the ONT 6mer model */
	for ( unsigned int i = 0; i < sequence.length() - 6; i++ ){

		sixMer = sequence.substr( i, 6 );

		if (i == BrdU_idx and isBrdU){

			nd.push_back( NormalDistribution( analogueModel[sixMer].first, analogueModel[sixMer].second ) );
		}
		else {

			nd.push_back( NormalDistribution( SixMer_model[sixMer].first, SixMer_model[sixMer].second ) );
		}
	}

	/*the first insertion state after start */
	State firstI = State( &ud, "-1_I", "", "", 1.0 );
	hmm.add_state( firstI );

	/*add states to the model, handle internal module transitions */
	for ( unsigned int i = 0; i < sequence.length() - 6; i++ ){

		loc = std::to_string( i );
		sixMer = sequence.substr( i, 6 );

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
	for ( unsigned int i = 0; i < sequence.length() - 7; i++ ){

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
	hmm.add_transition( states[0][sequence.length() - 7], hmm.end, 1.0 );
	hmm.add_transition( states[1][sequence.length() - 7], hmm.end, externalI2M1 );
	hmm.add_transition( states[2][sequence.length() - 7], hmm.end, externalM12M1 + externalM12D );

	hmm.finalise();

	return hmm.sequenceProbability( events );
}


int detect_main( int argc, char** argv ){

	Arguments trainArgs = parseDetectionArguments( argc, argv );

	/*import the analogue pore model that we specified on the command line */
	std::map< std::string, std::pair< double, double > > analogueModel =  import_poreModel( trainArgs.analogueModelFilename );

	/*get a filestream to the foh file - we'll load training data dynamically */
	std::ifstream fohFile( trainArgs.dataFilename );
	if ( not fohFile.is_open() ) throw IOerror( trainArgs.dataFilename );

	/*read the foh header - total count */
	std::string line;
	std::getline( fohFile, line );
	int detectionTotal = atoi(line.c_str());

	/*initialise progress */
	progressBar pb_align( detectionTotal );
	int prog = 0, offloadCount = 0, failed = 0;

	/*buffer */
	std::vector< read > buffer;

	/*open output file */
	std::ofstream outFile( trainArgs.outputFilename );
	if ( not outFile.is_open() ) throw IOerror( trainArgs.outputFilename );

	unsigned int windowLength = 20;

	/*align the training data */
	std::cout << "Detecting BrdU..." << std::endl;
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

		/*aligned pairs */
		std::getline( fohFile, line );
		std::istringstream sspairs( line );
		std::string qIndex, rIndex;
		while ( std::getline( sspairs, qIndex, ' ' ) ){

			std::getline( sspairs, rIndex, ' ' );

			currentRead.refToQuery[ atoi( rIndex.c_str() )] = atoi( qIndex.c_str() );
		}

		/*push it to the buffer or run Viterbi if the buffer is full */
		buffer.push_back( currentRead );
		if ( (buffer.size() < trainArgs.threads)  ) continue;

		/*HMM BrdU detection using the forward algorithm */
		#pragma omp parallel for schedule(dynamic) shared(failed,pb_align,buffer,trainArgs,detectionTotal,prog,SixMer_model,analogueModel) num_threads(trainArgs.threads)
		for ( auto r = buffer.begin(); r < buffer.end(); r++ ){

			/*normalise for shift and scale */
			eventDataForRead eventData = normaliseEvents( *r );

			//std::cout << eventData.qualityScore << std::endl;
			/*disregard this event if the quality score is too low */
			if ( fabs(eventData.qualityScore) > 2 ){

				failed++;
				prog++;
				continue;
			}

			/*get the subsequence of the reference this read mapped to and build an HMM from it */
			std::string refSeqMapped = (*r).mappedRefSubseq;
			
			/*push the filename for this read to the output */
			std::stringstream ss;
			ss << ">" << r -> filename << std::endl;

			int readHead = (r -> refToQuery)[0];

			/*exclude the starts and ends of the read, as the alignment tends to be worse there */
			for ( unsigned int i = windowLength; i < refSeqMapped.length() - 2*windowLength; i++ ){
			
				/*for each T we find in the read, calculate the log likelihood that it's an analogue */
				if ( analogueModel.count( refSeqMapped.substr(i, 6) ) > 0 ){

					int posOnQuery = (r -> refToQuery)[i];

					std::string readSnippet = refSeqMapped.substr(i - windowLength, 2*windowLength);
					std::vector< double > eventSnippet;

					/*get the events that correspond to the read snippet */
					for ( unsigned int j = readHead; j < (eventData.eventAlignment).size(); j++ ){

						/*move the readhead so that we can navigate through events efficiently */
						if ( (eventData.eventAlignment)[j].second == (r -> refToQuery)[i - windowLength] ) readHead = j;

						/*if an event has been aligned to a position in the window, add it */
						if ( (eventData.eventAlignment)[j].second >= (r -> refToQuery)[i - windowLength] and (eventData.eventAlignment)[j].second <= (r -> refToQuery)[i + windowLength] ){

							eventSnippet.push_back( (eventData.normalisedEvents)[j] );
						}

						/*stop once we get to the end of the window */
						if ( (eventData.eventAlignment)[j].second > (r -> refToQuery)[i + windowLength] ) break;
					}
					double logProbThymidine = seqProbability(readSnippet, eventSnippet, analogueModel, windowLength, false);
					double logProbAnalogue = seqProbability(readSnippet, eventSnippet, analogueModel, windowLength, true);
					double logLikelihoodRatio = logProbAnalogue - logProbThymidine;

					ss << i << "\t" << logLikelihoodRatio << "\t" << refSeqMapped.substr(i, 6) << "\t" << (r->basecall).substr(posOnQuery, 6) << std::endl;
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
		if ( fohFile.eof() ) break;
	}
	displayProgress( detectionTotal, detectionTotal );
	std::cout << "\nDone." << std::endl;

	outFile.close();
	fohFile.close();

	return 0;
}
