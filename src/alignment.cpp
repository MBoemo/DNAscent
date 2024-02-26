//----------------------------------------------------------
// Copyright 2019-2020 University of Oxford
// Written by Michael A. Boemo (mb915@cam.ac.uk)
// This software is licensed under GPL-3.0.  You should have
// received a copy of the license with this software.  If
// not, please Email the author.
//----------------------------------------------------------

//#define TEST_DETECT 1
//#define TEST_VITERBI 1

#include <fstream>
#include "detect.h"
#include <math.h>
#include <utility>
#include <stdlib.h>
#include <limits>
#include "common.h"
#include "event_handling.h"
#include "../fast5/include/fast5.hpp"
#include "alignment.h"
#include "error_handling.h"
#include "probability.h"
#include "htsInterface.h"


static const char *help=
"align: DNAscent executable that generates a BrdU-aware event alignment.\n"
"To run DNAscent align, do:\n"
"   DNAscent align -b /path/to/alignment.bam -r /path/to/reference.fasta -i /path/to/index.dnascent -o /path/to/output.detect\n"
"Required arguments are:\n"
"  -b,--bam                  path to alignment BAM file,\n"
"  -r,--reference            path to genome reference in fasta format,\n"
"  -i,--index                path to DNAscent index,\n"
"  -o,--output               path to output file that will be generated.\n"
"Optional arguments are:\n"
"  -t,--threads              number of threads (default is 1 thread),\n"
"  -m,--maxReads             maximum number of reads to consider,\n"
"  -q,--quality              minimum mapping quality (default is 20),\n"
"  -l,--length               minimum read length in bp (default is 100).\n"
"DNAscent is under active development by the Boemo Group, Department of Pathology, University of Cambridge (https://www.boemogroup.org/).\n"
"Please submit bug reports to GitHub Issues (https://github.com/MBoemo/DNAscent/issues).";

struct Arguments {
	std::string bamFilename;
	std::string referenceFilename;
	std::string outputFilename;
	std::string indexFilename;
	bool capReads;
	int minQ, maxReads;
	int minL;
	unsigned int threads;
};

Arguments parseAlignArguments( int argc, char** argv ){

	if( argc < 2 ){

		std::cout << "Exiting with error.  Insufficient arguments passed to DNAscent detect." << std::endl << help << std::endl;
		exit(EXIT_FAILURE);
	}

	if ( std::string( argv[ 1 ] ) == "-h" or std::string( argv[ 1 ] ) == "--help" ){

		std::cout << help << std::endl;
		exit(EXIT_SUCCESS);
	}
	else if( argc < 4 ){

		std::cout << "Exiting with error.  Insufficient arguments passed to DNAscent detect." << std::endl;
		exit(EXIT_FAILURE);
	}

	Arguments args;

	/*defaults - we'll override these if the option was specified by the user */
	args.threads = 1;
	args.minQ = 20;
	args.minL = 100;
	args.capReads = false;
	args.maxReads = 0;

	/*parse the command line arguments */

	for ( int i = 1; i < argc; ){

		std::string flag( argv[ i ] );

		if ( flag == "-b" or flag == "--bam" ){

			std::string strArg( argv[ i + 1 ] );
			args.bamFilename = strArg;
			i+=2;
		}
		else if ( flag == "-r" or flag == "--reference" ){

			std::string strArg( argv[ i + 1 ] );
			args.referenceFilename = strArg;
			i+=2;
		}
		else if ( flag == "-t" or flag == "--threads" ){

			std::string strArg( argv[ i + 1 ] );
			args.threads = std::stoi( strArg.c_str() );
			i+=2;
		}
		else if ( flag == "-q" or flag == "--quality" ){

			std::string strArg( argv[ i + 1 ] );
			args.minQ = std::stoi( strArg.c_str() );
			i+=2;
		}
		else if ( flag == "-l" or flag == "--length" ){

			std::string strArg( argv[ i + 1 ] );
			args.minL = std::stoi( strArg.c_str() );
			i+=2;
		}
		else if ( flag == "-i" or flag == "--index" ){

			std::string strArg( argv[ i + 1 ] );
			args.indexFilename = strArg;
			i+=2;
		}
		else if ( flag == "-o" or flag == "--output" ){

			std::string strArg( argv[ i + 1 ] );
			args.outputFilename = strArg;
			i+=2;
		}
		else if ( flag == "-m" or flag == "--maxReads" ){

			std::string strArg( argv[ i + 1 ] );
			args.capReads = true;
			args.maxReads = std::stoi( strArg.c_str() );
			i+=2;
		}
		else throw InvalidOption( flag );
	}
	if (args.outputFilename == args.indexFilename or args.outputFilename == args.referenceFilename or args.outputFilename == args.bamFilename) throw OverwriteFailure();

	return args;
}


inline double lnVecMax(std::vector<double> v){

	double maxVal = v[0];
	for (size_t i = 1; i < v.size(); i++){
		if (lnGreaterThan( v[i], maxVal )){
			maxVal = v[i];
		}
	}
	return maxVal;
}


inline int lnArgMax(std::vector<double> v){

	double maxVal = v[0];
	int maxarg = 0;

	for (size_t i = 1; i < v.size(); i++){
		if (lnGreaterThan( v[i], maxVal )){
			maxVal = v[i];
			maxarg = i;
		}
	}
	return maxarg;
}


std::pair< double, std::vector< std::string > > builtinViterbi( std::vector <double> &observations,
				std::string &sequence,
				PoreParameters scalings,
				bool flip){
	
	//HMM transition probabilities	
	double externalD2D = eln(Pore_Substrate_Config.HMM_config.externalD2D);
	double externalD2M1 = eln(Pore_Substrate_Config.HMM_config.externalD2M);
	double externalI2M1 = eln(Pore_Substrate_Config.HMM_config.externalI2M);
	double externalM12D = eln(Pore_Substrate_Config.HMM_config.externalM2D);
	double internalM12I = eln(Pore_Substrate_Config.HMM_config.internalM2I);
	double internalI2I = eln(Pore_Substrate_Config.HMM_config.internalI2I);

	//transition probabilities that change on a per-read basis
	double internalM12M1 = eln(1. - (1./scalings.eventsPerBase));
	double externalM12M1 = eln(1.0 - externalM12D - internalM12I - internalM12M1);
	double externalM12M1orD = lnSum( externalM12M1, externalM12D );
	double externalOrInternalM12M1 = lnSum( externalM12M1, internalM12M1 );
	int maxindex;

	unsigned int k = Pore_Substrate_Config.kmer_len;

	size_t n_states = sequence.length() - k + 1;

	//pre-compute 6mer indices
	std::vector<unsigned int> kmerIndices;
	kmerIndices.reserve(n_states);
	for (size_t i = 0; i < n_states; i++){
		std::string kmer = sequence.substr(i, k);
		if (flip) std::reverse(kmer.begin(),kmer.end());
		kmerIndices.push_back(kmer2index(kmer,k));
	}

	std::vector< std::vector< size_t > > backtraceS( 3*n_states, std::vector< size_t >( observations.size() + 1 ) ); /*stores state indices for the Viterbi backtrace */
	std::vector< std::vector< size_t > > backtraceT( 3*n_states, std::vector< size_t >( observations.size() + 1 ) ); /*stores observation indices for the Viterbi backtrace */

	//reserve 0 for start
	ssize_t D_offset = 0;
	ssize_t M_offset = n_states;
	ssize_t I_offset = 2*n_states;

	std::vector< double > I_curr(n_states, NAN), D_curr(n_states, NAN), M_curr(n_states, NAN), I_prev(n_states, NAN), D_prev(n_states, NAN), M_prev(n_states, NAN);
	double start_curr = NAN, start_prev = 0.0;

	double matchProb, insProb;

	/*-----------INITIALISATION----------- */
	//transitions from the start state
	D_prev[0] = lnProd( start_prev, externalM12D );
	backtraceS[0 + D_offset][0] = -1;
	backtraceT[0 + D_offset][0] = 0;

	//account for transitions between deletion states before we emit the first observation
	for ( unsigned int i = 1; i < n_states; i++ ){

		D_prev[i] = D_prev[i-1] + externalD2D;
		backtraceS[i + D_offset][0] = i -1 + D_offset;
		backtraceT[i + D_offset][0] = 0;
	}


	/*-----------RECURSION----------- */
	/*complexity is O(T*N^2) where T is the number of observations and N is the number of states */
	double level_mu, level_sigma;
	for ( unsigned int t = 0; t < observations.size(); t++ ){

		std::fill( I_curr.begin(), I_curr.end(), NAN );
		std::fill( M_curr.begin(), M_curr.end(), NAN );
		std::fill( D_curr.begin(), D_curr.end(), NAN );

		std::pair<double,double> meanStd = Pore_Substrate_Config.pore_model[kmerIndices[0]];

		//uncomment if you scale model
		//level_mu = (scalings.shift + scalings.scale * meanStd.first);
		//level_sigma = meanStd.second;

		//uncomment if you scale events
		level_mu = meanStd.first;
		level_sigma = meanStd.second;

		matchProb = eln( normalPDF( level_mu, level_sigma, (observations[t] - scalings.shift)/scalings.scale ) );
		//matchProb = eln( cauchyPDF( level_mu, level_sigma, (observations[t] - scalings.shift)/scalings.scale ) );
		insProb = 0.0; //log(1) = 0; set probability equal to 1 and use transition probability as weighting

		//to the base 1 insertion
		I_curr[0] = lnVecMax({I_prev[0] + internalI2I + insProb ,
			                  M_prev[0] + internalM12I + insProb,
							  start_prev + internalM12I + insProb
		                     });
		maxindex = lnArgMax({I_prev[0] + internalI2I + insProb,
					         M_prev[0] + internalM12I + insProb,
							 start_prev + internalM12I + insProb*0.001
						     });
		switch(maxindex){
			case 0:
				backtraceS[0 + I_offset][t+1] = 0 + I_offset;
				backtraceT[0 + I_offset][t+1] = t;
				break;
			case 1:
				backtraceS[0 + I_offset][t+1] = 0 + M_offset;
				backtraceT[0 + I_offset][t+1] = t;
				break;
			case 2:
				backtraceS[0 + I_offset][t+1] = -1;
				backtraceT[0 + I_offset][t+1] = t;
				break;
			default:
				std::cout << "problem" << std::endl;
				exit(EXIT_FAILURE);
		}

		//to the base 1 match
		M_curr[0] = lnVecMax({M_prev[0] + internalM12M1 + matchProb,
							  start_prev + externalOrInternalM12M1 + matchProb
							 });
		maxindex = lnArgMax({M_prev[0] + internalM12M1 + matchProb,
							 start_prev + externalOrInternalM12M1 + matchProb
							});
		switch(maxindex){
			case 0:
				backtraceS[0 + M_offset][t+1] = 0 + M_offset;
				backtraceT[0 + M_offset][t+1] = t;
				break;
			case 1:
				backtraceS[0 + M_offset][t+1] = -1;
				backtraceT[0 + M_offset][t+1] = t;
				break;
			default:
				std::cout << "problem" << std::endl;
				exit(EXIT_FAILURE);
		}

		//to the base 1 deletion
		D_curr[0] = lnProd( NAN, externalM12D );  //start to D
		backtraceS[0 + D_offset][t+1] = -1;
		backtraceT[0 + D_offset][t+1] = t + 1;


		//the rest of the sequence
		for ( unsigned int i = 1; i < n_states; i++ ){

			insProb = 0.0; //log(1) = 0; set probability equal to 1 and use transition probability as weighting

			//get model parameters
			std::pair<double,double> meanStd = Pore_Substrate_Config.pore_model[kmerIndices[i]];

			//uncomment if you scale model
			//level_mu = scalings.shift + scalings.scale * meanStd.first;
			//level_sigma = meanStd.second;

			//uncomment if you scale events
			level_mu = meanStd.first;
			level_sigma = meanStd.second;

			matchProb = eln( normalPDF( level_mu, level_sigma, (observations[t] - scalings.shift)/scalings.scale ) );
			//matchProb = eln( cauchyPDF( level_mu, level_sigma, (observations[t] - scalings.shift)/scalings.scale ) );
			
			//to the insertion
			I_curr[i] = lnVecMax({I_prev[i] + internalI2I + insProb,
								  M_prev[i] + internalM12I + insProb
								 });
			maxindex = lnArgMax({I_prev[i] + internalI2I + insProb ,
							     M_prev[i] + internalM12I + insProb
								});
			switch(maxindex){
				case 0:
					backtraceS[i + I_offset][t+1] = i + I_offset;
					backtraceT[i + I_offset][t+1] = t;
					break;
				case 1:
					backtraceS[i + I_offset][t+1] = i + M_offset;
					backtraceT[i + I_offset][t+1] = t;
					break;
				default:
					std::cout << "problem" << std::endl;
					exit(EXIT_FAILURE);
			}

			//to the match
			M_curr[i] = lnVecMax({I_prev[i-1] + externalI2M1 + matchProb,
								   M_prev[i-1] + externalM12M1 + matchProb ,
								   M_prev[i] + internalM12M1 + matchProb ,
								   D_prev[i-1] + externalD2M1 + matchProb
								   });
			maxindex = lnArgMax({I_prev[i-1] + externalI2M1 + matchProb,
							   M_prev[i-1] + externalM12M1 + matchProb,
							   M_prev[i] + internalM12M1 + matchProb,
							   D_prev[i-1] + externalD2M1 + matchProb
			});
			switch(maxindex){
				case 0:
					backtraceS[i + M_offset][t+1] = i - 1 + I_offset;
					backtraceT[i + M_offset][t+1] = t;
					break;
				case 1:
					backtraceS[i + M_offset][t+1] = i - 1 + M_offset;
					backtraceT[i + M_offset][t+1] = t;
					break;
				case 2:
					backtraceS[i + M_offset][t+1] = i + M_offset;
					backtraceT[i + M_offset][t+1] = t;
					break;
				case 3:
					backtraceS[i + M_offset][t+1] = i - 1 + D_offset;
					backtraceT[i + M_offset][t+1] = t;
					break;
				default:
					std::cout << "problem" << std::endl;
					exit(EXIT_FAILURE);
			}
		}

		for ( unsigned int i = 1; i < n_states; i++ ){

			//to the deletion
			D_curr[i] = lnVecMax({M_curr[i-1] + externalM12D,
				                  D_curr[i-1] + externalD2D
			                     });
			maxindex = lnArgMax({ M_curr[i-1] + externalM12D,
                               D_curr[i-1] + externalD2D
			                  });
			switch(maxindex){
				case 0:
					backtraceS[i + D_offset][t+1] = i - 1 + M_offset;
					backtraceT[i + D_offset][t+1] = t + 1;
					break;
				case 1:
					backtraceS[i + D_offset][t+1] = i - 1 + D_offset;
					backtraceT[i + D_offset][t+1] = t + 1;
					break;
				default:
					std::cout << "problem" << std::endl;
					exit(EXIT_FAILURE);
			}
		}

		I_prev = I_curr;
		M_prev = M_curr;
		D_prev = D_curr;
		start_prev = start_curr;
	}

	/*
	for (int i = 0; i < backtraceS.size(); i++){
		for (int j = 0; j < backtraceS[i].size(); j++){
			std::cout << backtraceS[i][j] << "\t";
		}
		std::cout << std::endl;
	}
	std::cout << "----------------------------------------------------------------------" << std::endl;
	*/

	/*-----------TERMINATION----------- */
	double viterbiScore = NAN;
	viterbiScore = lnVecMax( {D_curr.back() , // + eln( 1.0 ) which is 0 //D to end
							  M_curr.back() + externalM12M1orD,//M to end
							  I_curr.back() + externalI2M1 //I to end
							 });
	//std::cout << "Builtin Viterbi score: " << viterbiScore << std::endl;


	//figure out where to go from the end state
	maxindex = lnArgMax({ D_curr.back() , //+ eln( 1.0 ) which is 0
	                   M_curr.back() + externalM12M1orD ,
					   I_curr.back() + externalI2M1
	                   });

	ssize_t  traceback_new;
	ssize_t  traceback_old;
	ssize_t  traceback_t = observations.size();
	switch(maxindex){
		case 0:
			traceback_old = D_offset + n_states - 1;
			break;
		case 1:
			traceback_old = M_offset + n_states - 1;
			break;
		case 2:
			traceback_old = I_offset + n_states - 1;
			break;
		default:
			std::cout << "problem" << std::endl;
			exit(EXIT_FAILURE);
	}

#if TEST_VITERBI
std::cerr << "Starting traceback..." << std::endl;
std::cerr << "Number of events: " << observations.size() << std::endl;
#endif

	std::vector<std::string> stateIndices;
	stateIndices.reserve(observations.size());
	while (traceback_old != -1){

		traceback_new = backtraceS[ traceback_old ][ traceback_t ];
		traceback_t = backtraceT[ traceback_old ][ traceback_t ];

		if (traceback_old < M_offset){ //Del

			//std::cout << "D " << traceback_old << std::endl;
			stateIndices.push_back(std::to_string(traceback_old) + "_D");

		}
		else if (traceback_old < I_offset){ //M

			//std::cout << "M " << traceback_old - M_offset << std::endl;
			stateIndices.push_back(std::to_string(traceback_old- M_offset) + "_M");
		}
		else { //I

			//std::cout << "I " << traceback_old - I_offset << std::endl;
			stateIndices.push_back(std::to_string(traceback_old - I_offset) + "_I");

		}
		traceback_old = traceback_new;
	}
	std::reverse( stateIndices.begin(), stateIndices.end() );

#if TEST_VITERBI
std::cout << "Traceback terminated." << std::endl;
#endif

	return std::make_pair( viterbiScore,stateIndices);
}


bool referenceDefined(std::string &readSnippet){

	//make sure the read snippet is fully defined as A/T/G/C in reference
	unsigned int As = 0, Ts = 0, Cs = 0, Gs = 0;
	for ( std::string::iterator i = readSnippet.begin(); i < readSnippet.end(); i++ ){

		switch( *i ){
			case 'A' :
				As++;
				break;
			case 'T' :
				Ts++;
				break;
			case 'G' :
				Gs++;
				break;
			case 'C' :
				Cs++;
				break;
		}
	}
	if ( readSnippet.length() != (As + Ts + Gs + Cs) ){
		return false;
	}
	else return true;
}


std::shared_ptr<AlignedRead> eventalign( read &r, unsigned int totalWindowLength, std::map<unsigned int, std::pair<double,double>> &analogueCalls){

	//get the positions on the reference subsequence where we could attempt to make a call
	std::string strand;
	int readHead = 0;
	if ( r.isReverse ) strand = "rev";
	else strand = "fwd";

	unsigned int k = Pore_Substrate_Config.kmer_len;

	std::shared_ptr<AlignedRead> ar = std::make_shared<AlignedRead>(AlignedRead(r.readID, r.referenceMappedTo, strand, r.refStart, r.refEnd, (r.eventAlignment).size()));

	ar -> stdout += ">" + r.readID + " " + r.referenceMappedTo + " " + std::to_string(r.refStart) + " " + std::to_string(r.refEnd) + " " + strand + "\n";

	unsigned int posOnRef = 0;
	while ( posOnRef < r.referenceSeqMappedTo.size() ){

		//adjust so we can get the last bit of the read if it doesn't line up with the windows nicely
		unsigned int basesToEnd = r.referenceSeqMappedTo.size() - posOnRef ;
		unsigned int windowLength = std::min(basesToEnd, totalWindowLength);

		//find good breakpoints
		std::string break1, break2;
		if (basesToEnd > 1.5*totalWindowLength){

			std::string breakSnippet = (r.referenceSeqMappedTo).substr(posOnRef, 1.5*windowLength);

			bool isDefined = referenceDefined(breakSnippet);
			if (not isDefined){
				posOnRef += windowLength;
				continue;
			}

			for (unsigned int i = windowLength; i < 1.5*windowLength - k - 1; i++){

				std::string kmer = breakSnippet.substr(i,k);
				std::pair<double,double> meanStd = Pore_Substrate_Config.pore_model[kmer2index(kmer, k)];

				std::string kmer_back = breakSnippet.substr(i-1,k);
				std::pair<double,double> meanStd_back = Pore_Substrate_Config.pore_model[kmer2index(kmer_back, k)];

				std::string kmer_front = breakSnippet.substr(i+1,k);
				std::pair<double,double> meanStd_front = Pore_Substrate_Config.pore_model[kmer2index(kmer_front, k)];

				double gap1 = std::abs(meanStd.first - meanStd_front.first);
				double gap2 = std::abs(meanStd.first - meanStd_back.first);

				if (gap1 > 0.75 and gap2 > 0.75){
					break1 = breakSnippet.substr(i-1,k);
					break2 = breakSnippet.substr(i+1,k);
					windowLength = i + k;
					break;
				}
			}
		}

		std::string readSnippet = (r.referenceSeqMappedTo).substr(posOnRef, windowLength);

		bool isDefined = referenceDefined(readSnippet);
		if (not isDefined){

			posOnRef += windowLength;
			continue;
		}

		std::vector< double > eventSnippet_means;
		std::vector< event > eventSnippet;

		/*get the events that correspond to the read snippet */
		//ar -> stdout += "readHead at start: " + std::to_string(readHead) + "\n";
		bool firstMatch = true;
		int eventIdx = 0;
		for ( unsigned int j = readHead; j < (r.eventAlignment).size(); j++ ){

			/*if an event has been aligned to a position in the window, add it */
			if ( (r.refToQuery)[posOnRef] <= (r.eventAlignment)[j].second and (r.eventAlignment)[j].second < (r.refToQuery)[posOnRef + windowLength - k + 1] ){

				if (firstMatch){
					readHead = j;
					firstMatch = false;
				}

				std::vector<double> raw = (r.events)[(r.eventAlignment)[j].first].raw;
				eventSnippet_means.push_back((r.events)[(r.eventAlignment)[j].first].mean);
				eventSnippet.push_back((r.events)[(r.eventAlignment)[j].first]);
				eventIdx++;
			}

			/*stop once we get to the end of the window */
			if ( (r.eventAlignment)[j].second >= (r.refToQuery)[posOnRef + windowLength - k + 1] ) break;
		}
		
		//flag large insertions
		int querySpan = (r.refToQuery)[posOnRef + windowLength - k + 1] - (r.refToQuery)[posOnRef];
		assert(querySpan >= 0);
		int referenceSpan = windowLength - k + 1;
		int indelScore = querySpan - referenceSpan;

		//pass on this window if we have a deletion
		if ( eventSnippet_means.size() < 2){

			posOnRef += windowLength;
			continue;
		}

		//calculate where we are on the assembly - if we're a reverse complement, we're moving backwards down the reference genome
		int globalPosOnRef;
		if ( r.isReverse ) globalPosOnRef = r.refEnd - posOnRef - k;
		else globalPosOnRef = r.refStart + posOnRef;
		
		std::pair< double, std::vector<std::string> > builtinAlignment = builtinViterbi( eventSnippet_means, readSnippet, r.scalings, false);

		std::vector< std::string > stateLabels = builtinAlignment.second;
		size_t lastM_ev = 0;
		size_t lastM_ref = 0;

		size_t evIdx = 0;

		//grab the index of the last match so we don't print insertions where we shouldn't
		for (size_t i = 0; i < stateLabels.size(); i++){

			std::string label = stateLabels[i].substr(stateLabels[i].find('_')+1);
			int pos = std::stoi(stateLabels[i].substr(0,stateLabels[i].find('_')));

			if (label == "M"){
				lastM_ev = evIdx;
				lastM_ref = pos;
			}

			if (label != "D") evIdx++; //silent states don't emit an event
		}

		//do a second pass to print the alignment
		evIdx = 0;
		for (size_t i = 0; i < stateLabels.size(); i++){

			std::string label = stateLabels[i].substr(stateLabels[i].find('_')+1);
	        	int pos = std::stoi(stateLabels[i].substr(0,stateLabels[i].find('_')));

			if (label == "D") continue; //silent states don't emit an event

			std::string kmerStrand = (r.referenceSeqMappedTo).substr(posOnRef + pos, k);

			unsigned int evPos;
			std::string kmerRef;
			if (r.isReverse){
				evPos = globalPosOnRef - pos + k - 1;
				kmerRef = reverseComplement(kmerStrand);
			}
			else{
				evPos = globalPosOnRef + pos;
				kmerRef = kmerStrand;
			}

			if (label == "M"){
				std::pair<double,double> meanStd = Pore_Substrate_Config.pore_model[kmer2index(kmerStrand, k)];

				for (unsigned int idx_raw = 0; idx_raw < eventSnippet[evIdx].raw.size(); idx_raw++){
					double scaledEvent = (eventSnippet[evIdx].raw[idx_raw] - r.scalings.shift) / r.scalings.scale;
					if (analogueCalls.count(evPos) > 0){
						ar -> stdout += std::to_string(evPos) 
							      + "\t" + kmerRef 
							      + "\t" + std::to_string(scaledEvent) 
							      + "\t" + kmerStrand 
							      + "\t" + std::to_string(meanStd.first) 
							      + "\t" + std::to_string(analogueCalls.at(evPos).first) 
							      + "\t" + std::to_string(analogueCalls.at(evPos).second) 					          
							      + "\n";
					}
					else{
						ar -> stdout += std::to_string(evPos) 
							      + "\t" + kmerRef 
							      + "\t" + std::to_string(scaledEvent) 
							      + "\t" + kmerStrand 
							      + "\t" + std::to_string(meanStd.first) 
							      + "\n";
						ar -> addSignal(kmerStrand, evPos, scaledEvent, indelScore);
					}
				}

			}
			else if (label == "I" and evIdx < lastM_ev){ //don't print insertions after the last match because we're going to align these in the next segment
				for (unsigned int idx_raw = 0; idx_raw < eventSnippet[evIdx].raw.size(); idx_raw++){
					double scaledEvent = (eventSnippet[evIdx].raw[idx_raw] - r.scalings.shift) / r.scalings.scale;
					ar -> stdout += std::to_string(evPos) + "\t" + kmerRef + "\t" + std::to_string(scaledEvent) + "\t" + std::string(k, 'N') + "\t" + "0" + "\n";
				}
			}
			
			evIdx ++;
		}

		//TESTING - make sure nothing sketchy happens at the breakpoint
		//if (not found) out += "BREAKPOINT\n";
		//else ar.stdout += "BREAKPOINT PRIME " + break1 + " " + break2 + "\n";

		//go again starting at posOnRef + lastM_ref using events starting at readHead + lastM_ev
		//readHead += rawIdx2eventIdx[lastM_ev] + 1;
		readHead += lastM_ev + 1;
		posOnRef += lastM_ref + 1;
	}
	
	ar -> QCpassed = true;

	return ar;
}


int align_main( int argc, char** argv ){

	Arguments args = parseAlignArguments( argc, argv );

	//load DNAscent index
	std::map< std::string, std::string > readID2path;
	parseIndex( args.indexFilename, readID2path );

	//import fasta reference
	std::map< std::string, std::string > reference = import_reference_pfasta( args.referenceFilename );

	std::ofstream outFile( args.outputFilename );
	if ( not outFile.is_open() ) throw IOerror( args.outputFilename );

	htsFile* bam_fh;
	hts_idx_t* bam_idx;
	bam_hdr_t* bam_hdr;
	hts_itr_t* itr;

	//load the bam
	std::cout << "Opening bam file... ";
	bam_fh = sam_open((args.bamFilename).c_str(), "r");
	if (bam_fh == NULL) throw IOerror(args.bamFilename);

	//load the index
	bam_idx = sam_index_load(bam_fh, (args.bamFilename).c_str());
	if (bam_idx == NULL) throw IOerror("index for "+args.bamFilename);

	//load the header
	bam_hdr = sam_hdr_read(bam_fh);
	std::cout << "ok." << std::endl;

	/*initialise progress */
	int numOfRecords = 0, prog = 0, failed = 0;
	countRecords( bam_fh, bam_idx, bam_hdr, numOfRecords, args.minQ, args.minL );
	if (args.capReads){
		numOfRecords = std::min(numOfRecords,args.maxReads);
	}
	progressBar pb(numOfRecords,true);

	//build an iterator for all reads in the bam file
	const char *allReads = ".";
	itr = sam_itr_querys(bam_idx,bam_hdr,allReads);

	std::map<unsigned int, std::pair<double,double>> placeholder_analogueCalls;
	int result;
	int failedEvents = 0;
	unsigned int maxBufferSize;
	std::vector< bam1_t * > buffer;
	if ( args.threads <= 4 ) maxBufferSize = args.threads;
	else maxBufferSize = 4*(args.threads);

	do {
		//initialise the record and get the record from the file iterator
		bam1_t *record = bam_init1();
		result = sam_itr_next(bam_fh, itr, record);

		//add the record to the buffer if it passes the user's criteria, otherwise destroy it cleanly
		int mappingQual = record -> core.qual;
		int refStart,refEnd;
		getRefEnd(record,refStart,refEnd);
		int queryLen = record -> core.l_qseq;

		if ( mappingQual >= args.minQ and refEnd - refStart >= args.minL and queryLen != 0 ){

			buffer.push_back( record );
		}
		else{
			bam_destroy1(record);
		}

		/*if we've filled up the buffer with short reads, compute them in parallel */
		if (buffer.size() >= maxBufferSize or (buffer.size() > 0 and result == -1 ) ){

			#pragma omp parallel for schedule(dynamic) shared(buffer,Pore_Substrate_Config,args,prog,failed,placeholder_analogueCalls) num_threads(args.threads)
			for (unsigned int i = 0; i < buffer.size(); i++){

				read r;

				//get the read name (which will be the ONT readID from Albacore basecall)
				const char *queryName = bam_get_qname(buffer[i]);
				if (queryName == NULL) continue;
				std::string s_queryName(queryName);
				r.readID = s_queryName;

				//iterate on the cigar string to fill up the reference-to-query coordinate map
				parseCigar(buffer[i], r.refToQuery, r.queryToRef, r.refStart, r.refEnd);

				//get the name of the reference mapped to
				std::string mappedTo(bam_hdr -> target_name[buffer[i] -> core.tid]);
				r.referenceMappedTo = mappedTo;

				//open fast5 and normalise events to pA
				r.filename = readID2path[s_queryName];

				/*get the subsequence of the reference this read mapped to */
				r.referenceSeqMappedTo = reference.at(r.referenceMappedTo).substr(r.refStart, r.refEnd - r.refStart);

				//fetch the basecall from the bam file
				r.basecall = getQuerySequence(buffer[i]);

				//account for reverse complements
				if ( bam_is_rev(buffer[i]) ){

					r.basecall = reverseComplement( r.basecall );
					r.referenceSeqMappedTo = reverseComplement( r.referenceSeqMappedTo );
					r.isReverse = true;
				}

				bool useFitPoreModel = false;
				normaliseEvents(r, useFitPoreModel);

				//catch reads with rough event alignments that fail the QC
				if ( r.eventAlignment.size() == 0 ){

					failed++;
					prog++;
					continue;
				}

				std::shared_ptr<AlignedRead> ar = eventalign( r, Pore_Substrate_Config.windowLength_align, placeholder_analogueCalls);

				if (not ar -> QCpassed){
					failed++;
					prog++;
					continue;
				}

				#pragma omp critical
				{
					outFile << ar -> stdout;
					prog++;
					pb.displayProgress( prog, failed, failedEvents );
				}
			}
			for ( unsigned int i = 0; i < buffer.size(); i++ ) bam_destroy1(buffer[i]);
			buffer.clear();
		}
		pb.displayProgress( prog, failed, failedEvents );
		if (args.capReads and prog > args.maxReads){
			sam_itr_destroy(itr);
			return 0;
		}
	} while (result > 0);
	sam_itr_destroy(itr);
	std::cout << std::endl;
	return 0;
}
