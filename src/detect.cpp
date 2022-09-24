//----------------------------------------------------------
// Copyright 2019-2020 University of Oxford
// Written by Michael A. Boemo (mb915@cam.ac.uk)
// This software is licensed under GPL-3.0.  You should have
// received a copy of the license with this software.  If
// not, please Email the author.
//----------------------------------------------------------

//#define TEST_HMM 1
//#define TEST_LL 1
//#define TEST_ALIGNMENT 1
//#define TEST_METHYL 1

#include <fstream>
#include <algorithm>
#include <math.h>
#include <stdlib.h>
#include <limits>
#include "detect.h"
#include "common.h"
#include "event_handling.h"
#include "probability.h"
#include "../fast5/include/fast5.hpp"
#include "poreModels.h"
#include "../htslib/htslib/hts.h"
#include "../htslib/htslib/sam.h"
#include "../tensorflow/include/tensorflow/c/eager/c_api.h"
#include "htsInterface.h"
//#include "tensor.h"
//#include "alignment.h"
#include "error_handling.h"
#include <omp.h>


static const char *help=
"detect: DNAscent executable that detects BrdU in Oxford Nanopore reads.\n"
"To run DNAscent detect, do:\n"
"   DNAscent detect -b /path/to/alignment.bam -r /path/to/reference.fasta -i /path/to/index.dnascent -o /path/to/output.detect\n"
"Required arguments are:\n"
"  -b,--bam                  path to alignment BAM file,\n"
"  -r,--reference            path to genome reference in fasta format,\n"
"  -i,--index                path to DNAscent index,\n"
"  -o,--output               path to output file that will be generated.\n"
"Optional arguments are:\n"
"  -t,--threads              number of threads (default is 1 thread),\n"
"  --GPU                     use the GPU device indicated for prediction (default is CPU),\n"
"  -q,--quality              minimum mapping quality (default is 20),\n"
"  -l,--length               minimum read length in bp (default is 1000).\n"
"Written by Michael Boemo, Department of Pathology, University of Cambridge.\n"
"Please submit bug reports to GitHub Issues (https://github.com/MBoemo/DNAscent/issues).";

struct Arguments {
	std::string bamFilename;
	std::string referenceFilename;
	std::string outputFilename;
	std::string indexFilename;
	bool useGPU = false;
	unsigned char GPUdevice = '0';
	int minQ = 20;
	int minL = 1000;
	unsigned int threads = 1;
	double dilation = 1.0;
};

Arguments parseDetectArguments( int argc, char** argv ){

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
		else if ( flag == "--GPU" ){

			args.useGPU = true;
			std::string strArg( argv[ i + 1 ] );
			if (strArg.length() > 1) throw InvalidDevice(strArg);

			args.GPUdevice = *argv[ i + 1 ];

			i+=2;
		}
		else if ( flag == "--dilation" ){

			std::string strArg( argv[ i + 1 ] );
			args.dilation = std::stof( strArg.c_str() );
			i+=2;
		}
		else throw InvalidOption( flag );
	}
	if (args.outputFilename == args.indexFilename or args.outputFilename == args.referenceFilename or args.outputFilename == args.bamFilename) throw OverwriteFailure();

	return args;
}


double sequenceProbability( std::vector <double> &observations,
				std::string &sequence,
				size_t windowSize,
				bool useBrdU,
				PoreParameters scalings,
				size_t BrdUStart,
				size_t BrdUEnd ){
//covered in: tests/detect/hmm_forward

	//Initial transitions within modules (internal transitions)
	double internalM12I = eln(0.3475);
	double internalI2I = eln(0.5);
	double internalM12M1 = eln(0.4);

	//Initial transitions between modules (external transitions)
	double externalD2D = eln(0.3);
	double externalD2M1 = eln(0.7);
	double externalI2M1 = eln(0.5);
	double externalM12D = eln(0.0025);
	double externalM12M1 = eln(0.25);

	std::vector< double > I_curr(2*windowSize+1, NAN), D_curr(2*windowSize+1, NAN), M_curr(2*windowSize+1, NAN), I_prev(2*windowSize+1, NAN), D_prev(2*windowSize+1, NAN), M_prev(2*windowSize+1, NAN);
	double firstI_curr = NAN, firstI_prev = NAN;
	double start_curr = NAN, start_prev = 0.0;

	double matchProb, insProb;

	/*-----------INITIALISATION----------- */
	//transitions from the start state
	D_prev[0] = lnProd( start_prev, eln( 0.25 ) );

	//account for transitions between deletion states before we emit the first observation
	for ( unsigned int i = 1; i < D_prev.size(); i++ ){

		D_prev[i] = lnProd( D_prev[i-1], externalD2D );
	}


	/*-----------RECURSION----------- */
	/*complexity is O(T*N^2) where T is the number of observations and N is the number of states */
	double level_mu, level_sigma;
	for ( unsigned int t = 0; t < observations.size(); t++ ){

		std::fill( I_curr.begin(), I_curr.end(), NAN );
		std::fill( M_curr.begin(), M_curr.end(), NAN );
		std::fill( D_curr.begin(), D_curr.end(), NAN );
		firstI_curr = NAN;

		std::string sixMer = sequence.substr(0, 6);

		std::pair<double,double> meanStd = thymidineModel[sixMer2index(sixMer)];
		level_mu = scalings.shift + scalings.scale * meanStd.first;
		level_sigma = scalings.var * meanStd.second;

		//uncomment to scale events
		//level_mu = thymidineModel.at(sixMer).first;
		//level_sigma = scalings.var / scalings.scale * thymidineModel.at(sixMer).second;
		//observations[t] = (observations[t] - scalings.shift) / scalings.scale;

		matchProb = eln( normalPDF( level_mu, level_sigma, observations[t] ) );
		insProb = eln( uniformPDF( 0, 250, observations[t] ) );

		//first insertion
		firstI_curr = lnSum( firstI_curr, lnProd( lnProd( start_prev, eln( 0.25 ) ), insProb ) ); //start to first I
		firstI_curr = lnSum( firstI_curr, lnProd( lnProd( firstI_prev, eln( 0.25 ) ), insProb ) ); //first I to first I

		//to the base 1 insertion
		I_curr[0] = lnSum( I_curr[0], lnProd( lnProd( I_prev[0], internalI2I ), insProb ) );  //I to I
		I_curr[0] = lnSum( I_curr[0], lnProd( lnProd( M_prev[0], internalM12I ), insProb ) ); //M to I

		//to the base 1 match
		M_curr[0] = lnSum( M_curr[0], lnProd( lnProd( firstI_prev, eln( 0.5 ) ), matchProb ) ); //first I to first match
		M_curr[0] = lnSum( M_curr[0], lnProd( lnProd( M_prev[0], internalM12M1 ), matchProb ) );  //M to M
		M_curr[0] = lnSum( M_curr[0], lnProd( lnProd( start_prev, eln( 0.5 ) ), matchProb ) );  //start to M

		//to the base 1 deletion
		D_curr[0] = lnSum( D_curr[0], lnProd( NAN, eln( 0.25 ) ) );  //start to D
		D_curr[0] = lnSum( D_curr[0], lnProd( firstI_curr, eln( 0.25 ) ) ); //first I to first deletion

		//the rest of the sequence
		for ( unsigned int i = 1; i < I_curr.size(); i++ ){

			//get model parameters
			sixMer = sequence.substr(i, 6);
			std::pair<double,double> analogue_meanStd = analogueModel[sixMer2index(sixMer)];
			insProb = eln( uniformPDF( 0, 250, observations[t] ) );
			if ( useBrdU and BrdUStart <= i and i <= BrdUEnd and sixMer.find('T') != std::string::npos and analogue_meanStd.first != 0. ){

				level_mu = scalings.shift + scalings.scale * analogue_meanStd.first;
				level_sigma = scalings.var * analogue_meanStd.second;

				//uncomment if you scale events
				//level_mu = analogueModel.at(sixMer).first;
				//level_sigma = scalings.var / scalings.scale * analogueModel.at(sixMer).second;

				matchProb = eln( normalPDF( level_mu, level_sigma, observations[t] ) );
			}
			else{

				std::pair<double,double> meanStd = thymidineModel[sixMer2index(sixMer)];
				level_mu = scalings.shift + scalings.scale * meanStd.first;
				level_sigma = scalings.var * meanStd.second;

				//uncomment if you scale events
				//level_mu = thymidineModel.at(sixMer).first;
				//level_sigma = scalings.var / scalings.scale * thymidineModel.at(sixMer).second;

				matchProb = eln( normalPDF( level_mu, level_sigma, observations[t] ) );
			}

			//to the insertion
			I_curr[i] = lnSum( I_curr[i], lnProd( lnProd( I_prev[i], internalI2I ), insProb ) );  //I to I
			I_curr[i] = lnSum( I_curr[i], lnProd( lnProd( M_prev[i], internalM12I ), insProb ) ); //M to I

			//to the match
			M_curr[i] = lnSum( M_curr[i], lnProd( lnProd( I_prev[i-1], externalI2M1 ), matchProb ) );  //external I to M
			M_curr[i] = lnSum( M_curr[i], lnProd( lnProd( M_prev[i-1], externalM12M1 ), matchProb ) );  //external M to M
			M_curr[i] = lnSum( M_curr[i], lnProd( lnProd( M_prev[i], internalM12M1 ), matchProb ) );  //interal M to M
			M_curr[i] = lnSum( M_curr[i], lnProd( lnProd( D_prev[i-1], externalD2M1 ), matchProb ) );  //external D to M
		}

		for ( unsigned int i = 1; i < I_curr.size(); i++ ){

			//to the deletion
			D_curr[i] = lnSum( D_curr[i], lnProd( M_curr[i-1], externalM12D ) );  //external M to D
			D_curr[i] = lnSum( D_curr[i], lnProd( D_curr[i-1], externalD2D ) );  //external D to D
		}

		I_prev = I_curr;
		M_prev = M_curr;
		D_prev = D_curr;
		firstI_prev = firstI_curr;
		start_prev = start_curr;
	}


	/*-----------TERMINATION----------- */
	double forwardProb = NAN;

	forwardProb = lnSum( forwardProb, lnProd( D_curr.back(), eln( 1.0 ) ) ); //D to end
	forwardProb = lnSum( forwardProb, lnProd( M_curr.back(), lnSum(externalM12M1, externalM12D) ) ); //M to end
	forwardProb = lnSum( forwardProb, lnProd( I_curr.back(), externalI2M1 ) ); //I to end

#if TEST_HMM
std::cerr << "<-------------------" << std::endl;
std::cerr << useBrdU << std::endl;
std::cerr << scalings.shift << " " << scalings.scale << " " << scalings.var << std::endl;
std::cerr << sequence << std::endl;
for (auto ob = observations.begin(); ob < observations.end(); ob++){
	std::cerr << *ob << " ";
}
std::cerr << std::endl;
std::cerr << forwardProb << std::endl;
#endif

	return forwardProb;
}


void parseIndex( std::string indexFilename, std::map< std::string, std::string > &readID2path, bool &bulk ){

	std::cout << "Loading DNAscent index... ";
	std::ifstream indexFile( indexFilename );
	if ( not indexFile.is_open() ) throw IOerror( indexFilename );
	std::string line;

	//get whether this is bulk fast5 or individual fast5 from the index
	std::getline( indexFile, line);
	if (line == "#bulk") bulk = true;
	else if (line == "#individual") bulk = false;
	else throw IndexFormatting();

	//get the readID to path map
	while ( std::getline( indexFile, line) ){

		std::string readID = line.substr(0, line.find('\t'));
		std::string path = line.substr(line.find('\t')+1);
		readID2path[readID] = path;
	}
	std::cout << "ok." << std::endl;
}


std::vector< unsigned int > getPOIs( std::string &refSeq, int windowLength ){

	std::vector< unsigned int > POIs;

	for ( unsigned int i = 2*windowLength; i < refSeq.length() - 2*windowLength; i++ ){

		if (refSeq.substr(i,1) == "T") POIs.push_back(i);
	}
	return POIs;
}


std::string llAcrossRead( read &r,
                          unsigned int windowLength,
                          int &failedEvents,
                          bool methylAware ){

	std::string out;
	//get the positions on the reference subsequence where we could attempt to make a call
	std::vector< unsigned int > POIs = getPOIs( r.referenceSeqMappedTo, windowLength );
	std::string strand;
	unsigned int readHead = 0;
	if ( r.isReverse ){

		strand = "rev";
		readHead = (r.eventAlignment).size() - 1;
		std::reverse( POIs.begin(), POIs.end() );
	}
	else{

		strand = "fwd";
		readHead = 0;
	}

	out += ">" + r.readID + " " + r.referenceMappedTo + " " + std::to_string(r.refStart) + " " + std::to_string(r.refEnd) + " " + strand + "\n";

	for ( unsigned int i = 0; i < POIs.size(); i++ ){

		int posOnRef = POIs[i];
		int posOnQuery = (r.refToQuery).at(posOnRef);

		//sequence needs to be 6 bases longer than the span of events we catch
		//so sequence goes from posOnRef - windowLength to posOnRef + windowLength + 6
		//event span goes from posOnRef - windowLength to posOnRef + windowLength

		std::string readSnippet = (r.referenceSeqMappedTo).substr(posOnRef - windowLength, 2*windowLength+6);

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
		if ( readSnippet.length() != (As + Ts + Gs + Cs) ) continue;

		std::vector< double > eventSnippet;

		//catch spans with lots of insertions or deletions (this QC was set using results of tests/detect/hmm_falsePositives)
		unsigned int spanOnQuery = (r.refToQuery)[posOnRef + windowLength+6] - (r.refToQuery)[posOnRef - windowLength];
		if ( spanOnQuery > 3.5*windowLength or spanOnQuery < 2*windowLength ) continue;

		/*get the events that correspond to the read snippet */
		bool first = true;
		if ( r.isReverse ){

			for ( unsigned int j = readHead; j >= 0; j-- ){

				/*if an event has been aligned to a position in the window, add it */
				if ( (r.eventAlignment)[j].second >= (r.refToQuery)[posOnRef - windowLength] and (r.eventAlignment)[j].second < (r.refToQuery)[posOnRef + windowLength] ){

					if (first){
						readHead = j;
						first = false;
						//std::cout << "READHEAD:" << j << " " << readHead << std::endl;
					}

					double ev = (r.normalisedEvents)[(r.eventAlignment)[j].first];
					if (ev > 1.0 and ev < 250.0){
						eventSnippet.push_back( ev );
					}
					else{

						failedEvents++;
					}
				}

				/*stop once we get to the end of the window */
				if ( (r.eventAlignment)[j].second < (r.refToQuery)[posOnRef - windowLength] ){

					std::reverse(eventSnippet.begin(), eventSnippet.end());
					break;
				}
			}
		}
		else{
			for ( unsigned int j = readHead; j < (r.eventAlignment).size(); j++ ){

				/*if an event has been aligned to a position in the window, add it */
				if ( (r.eventAlignment)[j].second >= (r.refToQuery)[posOnRef - windowLength] and (r.eventAlignment)[j].second < (r.refToQuery)[posOnRef + windowLength] ){

					if (first){
						readHead = j;
						first = false;
						//std::cout << "READHEAD:" << j << " " << readHead << std::endl;
					}

					double ev = (r.normalisedEvents)[(r.eventAlignment)[j].first];
					if (ev > 1.0 and ev < 250.0){
						eventSnippet.push_back( ev );
					}
					else{

						failedEvents++;
					}
				}

				/*stop once we get to the end of the window */
				if ( (r.eventAlignment)[j].second >= (r.refToQuery)[posOnRef + windowLength] ) break;
			}
		}

		//catch abnormally few or many events (this QC was set using results of tests/detect/hmm_falsePositives)
		if ( eventSnippet.size() > 8*windowLength or eventSnippet.size() < 3.5*windowLength ) continue;

		/*
		TESTING - print out the read snippet, the ONT model, and the aligned events
		std::cout << readSnippet << std::endl;
		for ( int pos = 0; pos < readSnippet.length()-5; pos++ ){

			std::cout << readSnippet.substr(pos,6) << "\t" << thymidineModel.at( readSnippet.substr(pos,6) ).first << std::endl;
		}
		for ( auto ev = eventSnippet.begin(); ev < eventSnippet.end(); ev++){
			double scaledEv =  (*ev - r.scalings.shift) / r.scalings.scale;
			std::cout << scaledEv << std::endl;
		}
		*/

		//calculate where we are on the assembly - if we're a reverse complement, we're moving backwards down the reference genome
		int globalPosOnRef;
		std::string sixMerQuery = (r.basecall).substr(posOnQuery, 6);
		std::string sixMerRef = (r.referenceSeqMappedTo).substr(posOnRef, 6);
		if ( r.isReverse ){

			globalPosOnRef = r.refEnd - posOnRef - 1;
			sixMerQuery = reverseComplement( sixMerQuery );
			sixMerRef = reverseComplement( sixMerRef );
		}
		else{

			globalPosOnRef = r.refStart + posOnRef;
		}

		//make the BrdU call
		std::string sixOI = (r.referenceSeqMappedTo).substr(posOnRef,6);
		size_t BrdUStart = sixOI.find('T') + windowLength - 5;
		size_t BrdUEnd = windowLength;//sixOI.rfind('T') + windowLength;
		double logProbAnalogue = sequenceProbability( eventSnippet, readSnippet, windowLength, true, r.scalings, BrdUStart, BrdUEnd );
		double logProbThymidine = sequenceProbability( eventSnippet, readSnippet, windowLength, false, r.scalings, 0, 0 );
		double logLikelihoodRatio = logProbAnalogue - logProbThymidine;

#if TEST_LL
double runningKL = 0.0;
for (unsigned int s = 0; s < readSnippet.length() - 6; s++){
	std::string sixMer = readSnippet.substr(s,6);
	if ( BrdUStart <= s and s <= BrdUEnd and sixMer.find('T') != std::string::npos and analogueModel.count(sixMer) > 0 ){
		runningKL += KLdivergence( thymidineModel.at(sixMer).first, thymidineModel.at(sixMer).second, analogueModel.at(sixMer).first, analogueModel.at(sixMer).second );
	}
}
std::cerr << "<-------------------" << std::endl;
std::cerr << runningKL << std::endl;
std::cerr << spanOnQuery << std::endl;
std::cerr << readSnippet << std::endl;
for (auto ob = eventSnippet.begin(); ob < eventSnippet.end(); ob++){
	std::cerr << *ob << " ";
}
std::cerr << std::endl;
std::cerr << logLikelihoodRatio << std::endl;
#endif

		out += std::to_string(globalPosOnRef) + "\t" + std::to_string(logLikelihoodRatio) + "\t" + sixMerRef + "\t" + sixMerQuery + "\n";
	}
	return out;
}


std::map<unsigned int, double> llAcrossRead_forTraining( read &r, unsigned int windowLength){

	std::map<unsigned int, double> refPos2likelihood;

	//get the positions on the reference subsequence where we could attempt to make a call
	std::vector< unsigned int > POIs = getPOIs( r.referenceSeqMappedTo, windowLength );
	std::string strand;
	unsigned int readHead = 0;
	if ( r.isReverse ){

		strand = "rev";
		readHead = (r.eventAlignment).size() - 1;
		std::reverse( POIs.begin(), POIs.end() );
	}
	else{

		strand = "fwd";
		readHead = 0;
	}

	for ( unsigned int i = 0; i < POIs.size(); i++ ){

		int posOnRef = POIs[i];
		int posOnQuery = (r.refToQuery).at(posOnRef);

		//sequence needs to be 6 bases longer than the span of events we catch
		//so sequence goes from posOnRef - windowLength to posOnRef + windowLength + 6
		//event span goes from posOnRef - windowLength to posOnRef + windowLength

		std::string readSnippet = (r.referenceSeqMappedTo).substr(posOnRef - windowLength, 2*windowLength+6);

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
		if ( readSnippet.length() != (As + Ts + Gs + Cs) ) continue;

		//calculate where we are on the assembly - if we're a reverse complement, we're moving backwards down the reference genome
		int globalPosOnRef;
		std::string sixMerQuery = (r.basecall).substr(posOnQuery, 6);
		std::string sixMerRef = (r.referenceSeqMappedTo).substr(posOnRef, 6);
		if ( r.isReverse ){

			globalPosOnRef = r.refEnd - posOnRef - 6;
			sixMerQuery = reverseComplement( sixMerQuery );
			sixMerRef = reverseComplement( sixMerRef );
		}
		else{

			globalPosOnRef = r.refStart + posOnRef;
		}


		std::vector< double > eventSnippet;

		//catch spans with lots of insertions or deletions (this QC was set using results of tests/detect/hmm_falsePositives)
		unsigned int spanOnQuery = (r.refToQuery)[posOnRef + windowLength+6] - (r.refToQuery)[posOnRef - windowLength];
		if ( spanOnQuery > 3.5*windowLength or spanOnQuery < 2*windowLength ){
			refPos2likelihood[globalPosOnRef] = -20000; //tag an abort based on query span
			continue;
		}

		/*get the events that correspond to the read snippet */
		bool first = true;
		if ( r.isReverse ){

			for ( unsigned int j = readHead; j >= 0; j-- ){

				/*if an event has been aligned to a position in the window, add it */
				if ( (r.eventAlignment)[j].second >= (r.refToQuery)[posOnRef - windowLength] and (r.eventAlignment)[j].second < (r.refToQuery)[posOnRef + windowLength] ){

					if (first){
						readHead = j;
						first = false;
						//std::cout << "READHEAD:" << j << " " << readHead << std::endl;
					}

					double ev = (r.normalisedEvents)[(r.eventAlignment)[j].first];
					if (ev > 1.0 and ev < 250.0){
						eventSnippet.push_back( ev );
					}
				}

				/*stop once we get to the end of the window */
				if ( (r.eventAlignment)[j].second < (r.refToQuery)[posOnRef - windowLength] ){

					std::reverse(eventSnippet.begin(), eventSnippet.end());
					break;
				}
			}
		}
		else{
			for ( unsigned int j = readHead; j < (r.eventAlignment).size(); j++ ){

				/*if an event has been aligned to a position in the window, add it */
				if ( (r.eventAlignment)[j].second >= (r.refToQuery)[posOnRef - windowLength] and (r.eventAlignment)[j].second < (r.refToQuery)[posOnRef + windowLength] ){

					if (first){
						readHead = j;
						first = false;
						//std::cout << "READHEAD:" << j << " " << readHead << std::endl;
					}

					double ev = (r.normalisedEvents)[(r.eventAlignment)[j].first];
					if (ev > 1.0 and ev < 250.0){
						eventSnippet.push_back( ev );
					}
				}

				/*stop once we get to the end of the window */
				if ( (r.eventAlignment)[j].second >= (r.refToQuery)[posOnRef + windowLength] ) break;
			}
		}

		//make the BrdU call
		std::string sixOI = (r.referenceSeqMappedTo).substr(posOnRef,6);
		size_t BrdUStart = sixOI.find('T') + windowLength - 5;
		size_t BrdUEnd = windowLength;//sixOI.rfind('T') + windowLength;
		double logProbAnalogue = sequenceProbability( eventSnippet, readSnippet, windowLength, true, r.scalings, BrdUStart, BrdUEnd );
		double logProbThymidine = sequenceProbability( eventSnippet, readSnippet, windowLength, false, r.scalings, 0, 0 );
		double logLikelihoodRatio = logProbAnalogue - logProbThymidine;

		//catch abnormally few or many events (this QC was set using results of tests/detect/hmm_falsePositives)
		if ( eventSnippet.size() < 3.5*windowLength ){

			refPos2likelihood[globalPosOnRef] = -10000;//tag an abort based on number of events
			continue;
		}

#if TEST_LL
double runningKL = 0.0;
for (unsigned int s = 0; s < readSnippet.length() - 6; s++){
	std::string sixMer = readSnippet.substr(s,6);
	if ( BrdUStart <= s and s <= BrdUEnd and sixMer.find('T') != std::string::npos and analogueModel.count(sixMer) > 0 ){
		runningKL += KLdivergence( thymidineModel.at(sixMer).first, thymidineModel.at(sixMer).second, analogueModel.at(sixMer).first, analogueModel.at(sixMer).second );
	}
}
std::cerr << "<-------------------" << std::endl;
std::cerr << runningKL << std::endl;
std::cerr << spanOnQuery << std::endl;
std::cerr << readSnippet << std::endl;
for (auto ob = eventSnippet.begin(); ob < eventSnippet.end(); ob++){
	std::cerr << *ob << " ";
}
std::cerr << std::endl;
std::cerr << logLikelihoodRatio << std::endl;
#endif

		refPos2likelihood[globalPosOnRef] = logLikelihoodRatio;
	}
	return refPos2likelihood;
}


void read2tensor(std::shared_ptr<AlignedRead> r, const TensorShape &shape, TF_Tensor *t){

	std::vector<float> unformattedTensor = r -> makeTensor();

	size_t size = unformattedTensor.size();
	//put a check in here for size

	//auto output_array = std::make_unique<float[]>(size);
	std::cout << "mark1" << std::endl;
	float *output_array = (float *)malloc(size*sizeof(float));
	for(size_t i = 0; i < size; i++){
		output_array[i] = unformattedTensor[i];
	}
	std::cout << "mark2" << std::endl;

	t = TF_NewTensor(TF_FLOAT,
		shape.values,
		shape.dim,
		(void *)output_array,
		size*sizeof(float),
		cpp_array_deallocator<float>,
		nullptr);
	std::cout << "mark3" << std::endl;
	free(output_array);
}


std::string runCNN(std::shared_ptr<AlignedRead> r, std::shared_ptr<ModelSession> session){

	int NumInputs = 1;
	int NumOutputs = 1;

	std::vector<size_t> protoShape = r -> getShape();
	TensorShape input_shape={{1, (int64_t) protoShape[0], (int64_t) protoShape[1]}, 3};

	std::vector<float> unformattedTensor = r -> makeTensor();

	size_t size = unformattedTensor.size();
	assert(size > 0);

	float *tmp_array = (float *)malloc(size*sizeof(float));
	for(size_t i = 0; i < size; i++){
		tmp_array[i] = unformattedTensor[i];
	}

	TF_Tensor* InputValues = TF_NewTensor(TF_FLOAT,
		input_shape.values,
		input_shape.dim,
		(void *)tmp_array,
		size*sizeof(float),
		cpp_array_deallocator<float>,
		nullptr);

	TF_Tensor* OutputValues;

	//Run the Session
	CStatus status;
	TF_SessionRun(*(session->session.get()), NULL, &session->inputs, &InputValues, NumInputs, &session->outputs, &OutputValues, NumOutputs, NULL, 0, NULL, status.ptr);

	if(TF_GetCode(status.ptr) != TF_OK)
	{
		printf("%s",TF_Message(status.ptr));
	}

	if(TF_TensorType(OutputValues) != TF_FLOAT){
		std::cerr << "Error, unexpected output tensor type." << std::endl;
		exit (EXIT_FAILURE);
	}

	unsigned int outputFields = 3;

	//get positions on the read reference to write the output
	std::vector<unsigned int> positions = r -> getPositions();
	std::vector<int> alignmentQuality = r -> getAlignmentQuality();
	std::vector<std::string> sixMers = r -> getSixMers();

	size_t output_size = TF_TensorByteSize(OutputValues) / sizeof(float);
	assert(output_size == protoShape[0] * outputFields);
	float *output_array = (float *)TF_TensorData(OutputValues);

	//write the output
	unsigned int pos = 0;
	std::vector<std::string> lines;
	lines.reserve(positions.size());
	std::string thisPosition = std::to_string(positions[0]);
	std::string str_line, str_output;
	str_output += ">" + r -> getReadID() + " " + r -> getChromosome() + " " + std::to_string(r -> getMappingLower()) + " " + std::to_string(r -> getMappingUpper()) + " " + r -> getStrand() + "\n"; //header
	for(size_t i = 0; i < output_size; i++){
		if((i+1)%outputFields==0){

			//only output T positions
			if (sixMers[pos].substr(0,1) != "T"){
				pos++;
				continue;
			}

			str_line += thisPosition + "\t" + std::to_string(output_array[i])+ "\t" + std::to_string(output_array[i-1]);
			if (std::abs(alignmentQuality[pos]) > 100) str_line += "\t*";

			lines.push_back(str_line);
			str_line = "";
			pos++;
		}
		else{
			if (i != output_size-1) thisPosition = std::to_string(positions[pos]);
		}
	}

	TF_DeleteTensor(OutputValues);
	TF_DeleteTensor(InputValues);

	if (r -> getStrand() == "rev") std::reverse(lines.begin(),lines.end());

	for (auto s = lines.begin(); s < lines.end(); s++){
		str_output += *s + "\n";
	}

	return str_output;
}


std::map<unsigned int, std::pair<double,double>> runCNN_training(std::shared_ptr<AlignedRead> r, std::shared_ptr<ModelSession> session){

	int NumInputs = 1;
	int NumOutputs = 1;

	std::map<unsigned int, std::pair<double,double>> analogueCalls;

	std::vector<size_t> protoShape = r -> getShape();
	TensorShape input_shape={{1, (int64_t) protoShape[0], (int64_t) protoShape[1]}, 3};

	std::vector<float> unformattedTensor = r -> makeTensor();

	size_t size = unformattedTensor.size();
	assert(size > 0);

	float *tmp_array = (float *)malloc(size*sizeof(float));
	for(size_t i = 0; i < size; i++){
		tmp_array[i] = unformattedTensor[i];
	}

	TF_Tensor* InputValues = TF_NewTensor(TF_FLOAT,
		input_shape.values,
		input_shape.dim,
		(void *)tmp_array,
		size*sizeof(float),
		cpp_array_deallocator<float>,
		nullptr);

	TF_Tensor* OutputValues;

	//Run the Session
	CStatus status;
	TF_SessionRun(*(session->session.get()), NULL, &session->inputs, &InputValues, NumInputs, &session->outputs, &OutputValues, NumOutputs, NULL, 0, NULL, status.ptr);

	if(TF_GetCode(status.ptr) != TF_OK)
	{
		printf("%s",TF_Message(status.ptr));
	}

	if(TF_TensorType(OutputValues) != TF_FLOAT){
		std::cerr << "Error, unexpected output tensor type." << std::endl;
		exit (EXIT_FAILURE);
	}

	unsigned int outputFields = 3;

	//get positions on the read reference to write the output
	std::vector<unsigned int> positions = r -> getPositions();
	std::vector<std::string> sixMers = r -> getSixMers();

	size_t output_size = TF_TensorByteSize(OutputValues) / sizeof(float);
	assert(output_size == protoShape[0] * outputFields);
	float *output_array = (float *)TF_TensorData(OutputValues);

	//write the output
	unsigned int pos = 0;
	unsigned int thisPosition = positions[0];
	for(size_t i = 0; i < output_size; i++){
		if((i+1)%outputFields==0){

			//only output T positions
			if (sixMers[pos].substr(0,1) != "T"){
				pos++;
				continue;
			}

			analogueCalls[thisPosition] = std::make_pair(output_array[i],output_array[i-1]);
			pos++;
		}
		else{
			if (i != output_size-1) thisPosition = positions[pos];
		}
	}

	TF_DeleteTensor(OutputValues);
	TF_DeleteTensor(InputValues);

	return analogueCalls;
}


int detect_main( int argc, char** argv ){

	Arguments args = parseDetectArguments( argc, argv );
	bool bulkFast5;

	//load DNAscent index
	std::map< std::string, std::string > readID2path;
	parseIndex( args.indexFilename, readID2path, bulkFast5 );

	//get the neural network model path
	std::string pathExe = getExePath();
	std::string modelPath = pathExe + "dnn_models/detect_model_BrdUEdU/";
	std::string input_layer_name = "serving_default_input_1";

	std::shared_ptr<ModelSession> session;

	if (not args.useGPU){
		session = model_load_cpu(modelPath.c_str(), args.threads, input_layer_name.c_str());
	}
	else{
		session = model_load_gpu(modelPath.c_str(), args.GPUdevice, args.threads, input_layer_name.c_str());
	}

	//import fasta reference
	std::map< std::string, std::string > reference = import_reference_pfasta( args.referenceFilename );

	std::ofstream outFile( args.outputFilename );
	if ( not outFile.is_open() ) throw IOerror( args.outputFilename );

	//write the outfile header
	std::string outHeader = writeDetectHeader(args.bamFilename, args.referenceFilename, args.indexFilename, args.threads, false, args.minQ, args.minL, args.dilation, args.useGPU);
	outFile << outHeader;

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
	progressBar pb(numOfRecords,true);

	//build an iterator for all reads in the bam file
	const char *allReads = ".";
	itr = sam_itr_querys(bam_idx,bam_hdr,allReads);

	unsigned int windowLength_align = 50;

	int result;
	int failedEvents = 0;
	unsigned int maxBufferSize;
	std::vector< bam1_t * > buffer;
	maxBufferSize = 16*(args.threads);
	//maxBufferSize = args.threads;
	//if ( args.threads <= 4 ) maxBufferSize = args.threads;
	//else maxBufferSize = 4*(args.threads);

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

			#pragma omp parallel for schedule(dynamic) shared(buffer,windowLength_align,analogueModel,thymidineModel,args,prog,failed,session) num_threads(args.threads)
			for (unsigned int i = 0; i < buffer.size(); i++){

				read r;

				//get the read name (which will be the ONT readID from basecall)
				const char *queryName = bam_get_qname(buffer[i]);
				if (queryName == NULL) continue;
				std::string s_queryName(queryName);
				r.readID = s_queryName;

				//iterate on the cigar string to fill up the reference-to-query coordinate map
				parseCigar(buffer[i], r.refToQuery, r.refStart, r.refEnd);

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

				normaliseEvents(r, bulkFast5);

				//catch reads with rough event alignments that fail the QC
				if ( r.eventAlignment.size() == 0 ){

					failed++;
					prog++;
					continue;
				}

				std::pair<bool,std::shared_ptr<AlignedRead>> ar = eventalign_detect( r, windowLength_align, args.dilation );

				if (not ar.first){
					failed++;
					prog++;
					continue;
				}

				std::string readOut = runCNN(ar.second,session);

				prog++;
				pb.displayProgress( prog, failed, failedEvents );
				
				#pragma omp critical
				{
					outFile << readOut;
					pb.displayProgress( prog, failed, failedEvents );
				}
				

			}

			for ( unsigned int i = 0; i < buffer.size(); i++ ) bam_destroy1(buffer[i]);
			buffer.clear();
		}
		pb.displayProgress( prog, failed, failedEvents );
	} while (result > 0);
	sam_itr_destroy(itr);
	std::cout << std::endl;
	return 0;
}
