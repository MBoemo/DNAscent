//----------------------------------------------------------
// Copyright 2019 University of Oxford
// Written by Michael A. Boemo (michael.boemo@path.ox.ac.uk)
// This software is licensed under GPL-2.0.  You should have
// received a copy of the license with this software.  If
// not, please Email the author.
//----------------------------------------------------------

#include <fstream>
#include "detect.h"
#include <math.h>
#include <stdlib.h>
#include <limits>
#include "common.h"
#include "data_IO.h"
#include "event_handling.h"
#include "probability.h"
#include "../htslib/htslib/hts.h"
#include "../htslib/htslib/sam.h"
#include "../fast5/include/fast5.hpp"


static const char *help=
"detect: DNAscent executable that detects BrdU in Oxford Nanopore reads.\n"
"To run DNAscent detect, do:\n"
"  ./DNAscent detect [arguments]\n"
"Example:\n"
"  ./DNAscent detect -b /path/to/alignment.bam -r /path/to/reference.fasta -i /path/to/index.index -o /path/to/output.out -t 20\n"
"Required arguments are:\n"
"  -b,--bam                  path to alignment BAM file,\n"
"  -r,--reference            path to genome reference in fasta format,\n"
"  -i,--index                path to DNAscent index,\n"
"  -o,--output               path to output file that will be generated.\n"
"Optional arguments are:\n"
"  -t,--threads              number of threads (default is 1 thread)\n"
"  --divergence              minimum KL-divergence between BrdU 6mers to include and ONT pore model (default is 2.0),\n"
"  --noCpG                   exclude 6mers that contain a CpG (default is to include all 6mers),\n"
"  -q,--quality              minimum mapping quality (default is 20).\n"
"  -l,--length               minimum read length in bp (default is 100).\n";

struct Arguments {
	std::string bamFilename;
	std::string referenceFilename;
	std::string outputFilename;
	std::string indexFilename;
	bool excludeCpG;
	double divergence;
	int minQ;
	unsigned int minL;
	unsigned int threads;
};

extern std::map< std::string, std::pair< double, double > > SixMer_model;
extern std::map< std::string, std::pair< double, double > > BrdU_model_full;

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

	/*defaults - we'll override these if the option was specified by the user */
	args.threads = 1;
	args.minQ = 20;
	args.minL = 100;
	args.excludeCpG = false;
	args.divergence = 2.0;

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
		else if ( flag == "--divergence" ){

			std::string strArg( argv[ i + 1 ] );
			args.divergence = std::stof(strArg.c_str());
			i+=2;
		}
		else if ( flag == "--noCpG" ){

			args.excludeCpG = true;
			i+=1;
		}
		else throw InvalidOption( flag );
	}
	return args;
}


//Initial transitions within modules (internal transitions)
static double internalM12I = 0.001;
static double internalI2I = 0.001;
static double internalM12M1 = 0.4;

//Initial transitions between modules (external transitions)
static double externalD2D = 0.3;
static double externalD2M1 = 0.7;
static double externalI2M1 = 0.999;
static double externalM12D = 0.0025;
static double externalM12M1 = 0.5965;

double sequenceProbability( std::vector <double> &observations, std::string &sequence, size_t windowSize, bool useBrdU, std::map< std::string, std::pair< double, double > > &analogueModel, PoreParameters scalings ){

	std::vector< double > I_curr(2*windowSize-5, NAN), D_curr(2*windowSize-5, NAN), M_curr(2*windowSize-5, NAN), I_prev(2*windowSize-5, NAN), D_prev(2*windowSize-5, NAN), M_prev(2*windowSize-5, NAN);
	double firstI_curr = NAN, firstI_prev = NAN;
	double start_curr = NAN, start_prev = 0.0;

	double matchProb, insProb;

	/*-----------INITIALISATION----------- */
	//transitions from the start state
	D_prev[0] = lnProd( start_prev, eln( 0.25 ) );

	//account for transitions between deletion states before we emit the first observation
	for ( unsigned int i = 1; i < D_prev.size(); i++ ){

		D_prev[i] = lnProd( D_prev[i-1], eln ( externalD2D ) );
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

		level_mu = scalings.shift + scalings.scale * SixMer_model.at(sixMer).first;
		level_sigma = scalings.var * SixMer_model.at(sixMer).second;
		//level_sigma = SixMer_model.at(sixMer).second;

		//uncomment to scale events
		//level_mu = SixMer_model.at(sixMer).first;
		//level_sigma = scalings.var / scalings.scale * SixMer_model.at(sixMer).second;
		//observations[t] = (observations[t] - scalings.shift) / scalings.scale;

		matchProb = eln( normalPDF( level_mu, level_sigma, observations[t] ) );
		insProb = eln( uniformPDF( 0, 250, observations[t] ) );

		//first insertion
		firstI_curr = lnSum( firstI_curr, lnProd( lnProd( start_prev, eln( 0.25 ) ), insProb ) ); //start to first I
		firstI_curr = lnSum( firstI_curr, lnProd( lnProd( firstI_prev, eln( 0.25 ) ), insProb ) ); //first I to first I

		//to the base 1 insertion
		I_curr[0] = lnSum( I_curr[0], lnProd( lnProd( I_prev[0], eln( internalI2I ) ), insProb ) );  //I to I
		I_curr[0] = lnSum( I_curr[0], lnProd( lnProd( M_prev[0], eln( internalM12I ) ), insProb ) ); //M to I 

		//to the base 1 match
		M_curr[0] = lnSum( M_curr[0], lnProd( lnProd( firstI_prev, eln( 0.5 ) ), matchProb ) ); //first I to first match
		M_curr[0] = lnSum( M_curr[0], lnProd( lnProd( M_prev[0], eln( internalM12M1 ) ), matchProb ) );  //M to M
		M_curr[0] = lnSum( M_curr[0], lnProd( lnProd( start_prev, eln( 0.5 ) ), matchProb ) );  //start to M

		//to the base 1 deletion
		D_curr[0] = lnSum( D_curr[0], lnProd( NAN, eln( 0.25 ) ) );  //start to D
		D_curr[0] = lnSum( D_curr[0], lnProd( firstI_curr, eln( 0.25 ) ) ); //first I to first deletion

		//the rest of the sequence
		for ( unsigned int i = 1; i < I_curr.size(); i++ ){

			//get model parameters
			sixMer = sequence.substr(i, 6);
			insProb = eln( uniformPDF( 0, 250, observations[t] ) );
			if ( useBrdU and i == windowSize ){

				level_mu = scalings.shift + scalings.scale * analogueModel.at(sixMer).first;
				level_sigma = scalings.var * analogueModel.at(sixMer).second;
				//level_sigma = analogueModel.at(sixMer).second;

				//uncomment if you scale events
				//level_mu = analogueModel.at(sixMer).first;
				//level_sigma = scalings.var / scalings.scale * analogueModel.at(sixMer).second;

				matchProb = eln( normalPDF( level_mu, level_sigma, observations[t] ) );
			}
			else{

				level_mu = scalings.shift + scalings.scale * SixMer_model.at(sixMer).first;
				level_sigma = scalings.var * SixMer_model.at(sixMer).second;
				//level_sigma = SixMer_model.at(sixMer).second;

				//uncomment if you scale events				
				//level_mu = SixMer_model.at(sixMer).first;
				//level_sigma = scalings.var / scalings.scale * SixMer_model.at(sixMer).second;

				matchProb = eln( normalPDF( level_mu, level_sigma, observations[t] ) );
			}

			//to the insertion
			I_curr[i] = lnSum( I_curr[i], lnProd( lnProd( I_prev[i], eln( internalI2I ) ), insProb ) );  //I to I
			I_curr[i] = lnSum( I_curr[i], lnProd( lnProd( M_prev[i], eln( internalM12I ) ), insProb ) ); //M to I 

			//to the match
			M_curr[i] = lnSum( M_curr[i], lnProd( lnProd( I_prev[i-1], eln( externalI2M1 ) ), matchProb ) );  //external I to M
			M_curr[i] = lnSum( M_curr[i], lnProd( lnProd( M_prev[i-1], eln( externalM12M1 ) ), matchProb ) );  //external M to M
			M_curr[i] = lnSum( M_curr[i], lnProd( lnProd( M_prev[i], eln( internalM12M1 ) ), matchProb ) );  //interal M to M
			M_curr[i] = lnSum( M_curr[i], lnProd( lnProd( D_prev[i-1], eln( externalD2M1 ) ), matchProb ) );  //external D to M
		}

		for ( unsigned int i = 1; i < I_curr.size(); i++ ){

			//to the deletion
			D_curr[i] = lnSum( D_curr[i], lnProd( M_curr[i-1], eln( externalM12D ) ) );  //external M to D
			D_curr[i] = lnSum( D_curr[i], lnProd( D_curr[i-1], eln( externalD2D ) ) );  //external D to D
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
	forwardProb = lnSum( forwardProb, lnProd( M_curr.back(), eln( externalM12M1 + externalM12D ) ) ); //M to end
	forwardProb = lnSum( forwardProb, lnProd( I_curr.back(), eln( externalI2M1 ) ) ); //I to end

	return forwardProb;
}


std::string getQuerySequence( bam1_t *record ){ 
	
	std::string seq;
	uint8_t *a_seq = bam_get_seq(record);
	for ( int i = 0; i < record -> core.l_qseq; i++){
		int seqInBase = bam_seqi(a_seq,i);

		switch (seqInBase) {

			case 1: seq += "A"; break;
			case 2: seq += "C"; break;
			case 4: seq += "G"; break;
			case 8: seq += "T"; break;
			case 15: seq += "N"; break;
			default: throw ParsingError();
		}
	}
	return seq;
}


void parseCigar(bam1_t *record, std::map< unsigned int, unsigned int > &ref2query, int &refStart, int &refEnd ){

	//initialise reference and query coordinates for the first match
	refStart = record -> core.pos;
	int queryPosition = 0;
	int refPosition = 0;

	const uint32_t *cigar = bam_get_cigar(record);

	if ( bam_is_rev(record) ){

		for ( int i = record -> core.n_cigar - 1; i >= 0; i--){

			const int op = bam_cigar_op(cigar[i]); //cigar operation
			const int ol = bam_cigar_oplen(cigar[i]); //number of consecutive operations

			//for a match, advance both reference and query together
			if (op == BAM_CMATCH or op == BAM_CEQUAL or op == BAM_CDIFF){

				for ( int j = refPosition; j < refPosition + ol; j++ ){

					ref2query[j] = queryPosition;
					queryPosition++;
				}
				refPosition += ol;
			}
			//for a deletion, advance only the reference position
			else if (op == BAM_CDEL or op == BAM_CREF_SKIP){

				for ( int j = refPosition; j < refPosition + ol; j++ ){

					ref2query[j] = queryPosition;
				}
				refPosition += ol;
			}
			//for insertions or soft clipping, advance only the query position
			else if (op == BAM_CSOFT_CLIP or op == BAM_CINS){

				for ( int j = refPosition; j < refPosition + ol; j++ ){

					ref2query[j] = queryPosition;
					queryPosition++;
				}
			}
			//N.B. hard clipping advances neither refernce nor query, so ignore it
		}
	}
	else {

		for ( unsigned int i = 0; i < record -> core.n_cigar; ++i){

			const int op = bam_cigar_op(cigar[i]); //cigar operation
			const int ol = bam_cigar_oplen(cigar[i]); //number of consecutive operations

			//for a match, advance both reference and query together
			if (op == BAM_CMATCH or op == BAM_CEQUAL or op == BAM_CDIFF){

				for ( int j = refPosition; j < refPosition + ol; j++ ){

					ref2query[j] = queryPosition;
					queryPosition++;
				}
				refPosition += ol;
			}
			//for a deletion, advance only the reference position
			else if (op == BAM_CDEL or op == BAM_CREF_SKIP){

				for ( int j = refPosition; j < refPosition + ol; j++ ){

					ref2query[j] = queryPosition;
				}
				refPosition += ol;
			}
			//for insertions or soft clipping, advance only the query position
			else if (op == BAM_CSOFT_CLIP or op == BAM_CINS){

				for ( int j = refPosition; j < refPosition + ol; j++ ){

					ref2query[j] = queryPosition;
					queryPosition++;
				}
			}
			//N.B. hard clipping advances neither refernce nor query, so ignore it
		}
	}
	refEnd = refStart + refPosition;
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

void countRecords( htsFile *bam_fh, hts_idx_t *bam_idx, bam_hdr_t *bam_hdr, int &numOfRecords, int minQ, int minL ){

	hts_itr_t* itr = sam_itr_querys(bam_idx,bam_hdr,".");
	int result;

	do {
		bam1_t *record = bam_init1();
		result = sam_itr_next(bam_fh, itr, record);
		if ( (record -> core.qual >= minQ) and (record -> core.l_qseq >= minL) ) numOfRecords++;
		bam_destroy1(record);
	} while (result > 0);

	//cleanup
	sam_itr_destroy(itr);
}


std::vector< unsigned int > getPOIs( std::string &refSeq, std::map< std::string, std::pair< double, double > > &analogueModel, int windowLength ){

	std::vector< unsigned int > POIs;

	for ( unsigned int i = 2*windowLength; i < refSeq.length() - 2*windowLength; i++ ){

		if ( analogueModel.count( refSeq.substr(i, 6) ) > 0 ) POIs.push_back(i);
	}
	return POIs;
}


void llAcrossRead( read &r, unsigned int windowLength, std::map< std::string, std::pair< double, double > > &analogueModel, std::stringstream &ss ){

	//get the positions on the reference subsequence where we could attempt to make a call
	std::vector< unsigned int > POIs = getPOIs( r.referenceSeqMappedTo, analogueModel, windowLength );

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

	ss << ">" << r.readID << " " << r.referenceMappedTo << " " << r.refStart << " " << r.refEnd << " " << strand << std::endl;

	for ( unsigned int i = 0; i < POIs.size(); i++ ){

		int posOnRef = POIs[i];
		int posOnQuery = (r.refToQuery).at(posOnRef);

		std::string readSnippet = (r.referenceSeqMappedTo).substr(posOnRef - windowLength, 2*windowLength);

		//TESTING - print out the read snippet and the event and the ONT model
		/*
		std::cout << "ref start: " << r.refStart << std::endl;
		std::cout << "ref end: " << r.refEnd << std::endl;
		std::cout << "position on ref: " << posOnRef << std::endl;
		std::cout << "strand: " << strand << std::endl;
		extern std::map< std::string, std::pair< double, double > > SixMer_model;
		std::cout << readSnippet << std::endl;
		for ( int pos = 0; pos < readSnippet.length()-5; pos++ ){
		
			std::cout << readSnippet.substr(pos,6) << "\t" << SixMer_model.at( readSnippet.substr(pos,6) ).first << std::endl;
		}
		*/
		//END TESTING

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

		//catch spans with lots of insertions or deletions
		int spanOnQuery = (r.refToQuery)[posOnRef + windowLength] - (r.refToQuery)[posOnRef - windowLength];
		if ( spanOnQuery > 2.5*windowLength or spanOnQuery < 1.5*windowLength ) continue;

		/*get the events that correspond to the read snippet */
		bool first = true;
		if ( r.isReverse ){

			for ( unsigned int j = readHead; j >= 0; j-- ){

				/*if an event has been aligned to a position in the window, add it */
				if ( (r.eventAlignment)[j].second >= (r.refToQuery)[posOnRef - windowLength] and (r.eventAlignment)[j].second < (r.refToQuery)[posOnRef + windowLength - 5] ){

					if (first){
						readHead = j;
						first = false;
					}

					double ev = (r.normalisedEvents)[(r.eventAlignment)[j].first];
					if (ev > 0 and ev < 250) eventSnippet.push_back( ev );

					//TESTING - print the event snippet
					//std::cout << "snippet size: " << eventSnippet.size() << std::endl;
					//std::cout << j << " " << (r.eventAlignment)[j].first << " " << (r.eventAlignment)[j].second << std::endl;
					//std::cout << (ev - r.scalings.shift) / r.scalings.scale << std::endl;
					//END TESTING
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
				if ( (r.eventAlignment)[j].second >= (r.refToQuery)[posOnRef - windowLength] and (r.eventAlignment)[j].second < (r.refToQuery)[posOnRef + windowLength - 5] ){

					if (first){
						readHead = j;
						first = false;
					}

					double ev = (r.normalisedEvents)[(r.eventAlignment)[j].first];
					if (ev > 0 and ev < 250) eventSnippet.push_back( ev );

					//TESTING - print the event snippet
					//std::cout << "snippet size: " << eventSnippet.size() << std::endl;
					//std::cout << j << " " << (r.eventAlignment)[j].first << " " << (r.eventAlignment)[j].second << std::endl;
					//std::cout << (ev - r.scalings.shift) / r.scalings.scale << std::endl;
					//END TESTING
				}

				/*stop once we get to the end of the window */
				if ( (r.eventAlignment)[j].second >= (r.refToQuery)[posOnRef + windowLength - 5] ) break;
			}
		}

		//catch abnormally few or many events
		if ( eventSnippet.size() > 8*windowLength or eventSnippet.size() < windowLength ) continue;

		std::string sixOI = (r.referenceSeqMappedTo).substr(posOnRef,6);
		std::vector<double> BrdUscores;
		double logProbAnalogue = sequenceProbability( eventSnippet, readSnippet, windowLength, true, analogueModel, r.scalings );

		//std::cout << "log likelihood brdu: " << logProbAnalogue << std::endl;
		double logProbThymidine = sequenceProbability( eventSnippet, readSnippet, windowLength, false, analogueModel, r.scalings );
		//std::cout << "log likelihood thym: " << logProbThymidine << std::endl;
		double logLikelihoodRatio = logProbAnalogue - logProbThymidine;
		//std::cout << "log likelihood ratio:" << logLikelihoodRatio << std::endl;
		//std::cout << "----------------------------------------------" << std::endl;

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
		ss << globalPosOnRef << "\t" << logLikelihoodRatio << "\t" << sixMerRef << "\t" << sixMerQuery << std::endl;
	}
}


int detect_main( int argc, char** argv ){

	Arguments args = parseDetectArguments( argc, argv );
	bool bulkFast5;

	/*import the analogue pore model that we specified on the command line */
	std::map< std::string, std::pair< double, double > > analogueModel = buildAnalogueModel(args.divergence, args.excludeCpG);

	//load DNAscent index
	std::map< std::string, std::string > readID2path;
	parseIndex( args.indexFilename, readID2path, bulkFast5 );

	//import fasta reference
	std::map< std::string, std::string > reference = import_reference_pfasta( args.referenceFilename );

	std::ofstream outFile( args.outputFilename );
	if ( not outFile.is_open() ) throw IOerror( args.outputFilename );

	htsFile* bam_fh;
	hts_idx_t* bam_idx;
	bam_hdr_t* bam_hdr;
	hts_itr_t* itr;

	//load the bam
	bam_fh = sam_open((args.bamFilename).c_str(), "r");
	if (bam_fh == NULL) throw IOerror(args.bamFilename);

	//load the index
	bam_idx = sam_index_load(bam_fh, (args.bamFilename).c_str());
	if (bam_idx == NULL) throw IOerror("index for "+args.bamFilename);

	//load the header
	bam_hdr = sam_hdr_read(bam_fh);

	/*initialise progress */
	int numOfRecords = 0, prog = 0, failed = 0;
	countRecords( bam_fh, bam_idx, bam_hdr, numOfRecords, args.minQ, args.minL );
	progressBar pb(numOfRecords,true);

	//build an iterator for all reads in the bam file
	const char *allReads = ".";
	itr = sam_itr_querys(bam_idx,bam_hdr,allReads);

	unsigned int windowLength = 10;
	int result;
	unsigned int maxBufferSize;
	std::vector< bam1_t * > buffer;
	if ( args.threads <= 4 ) maxBufferSize = args.threads;
	else maxBufferSize = 4*(args.threads);

	do {

		//initialise the record and get the record from the file iterator
		bam1_t *record = bam_init1();
		result = sam_itr_next(bam_fh, itr, record);

		int mappingQual = record -> core.qual;

		//add the record to the buffer if it passes the user's criteria, otherwise destroy it cleanly
		if ( mappingQual >= args.minQ ){
			buffer.push_back( record );
		}
		else{
			bam_destroy1(record);
		}
		
		/*if we've filled up the buffer with short reads, compute them in parallel */
		if (buffer.size() >= maxBufferSize or (buffer.size() > 0 and result == -1 ) ){

			#pragma omp parallel for schedule(dynamic) shared(buffer,windowLength,analogueModel,args,prog,failed) num_threads(args.threads)
			for (unsigned int i = 0; i < buffer.size(); i++){

				read r; 

				//get the read name (which will be the ONT readID from Albacore basecall)
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
				try{

					if (bulkFast5) bulk_getEvents(r.filename, r.readID, r.raw);			
					else getEvents( r.filename, r.raw );
				}
				catch ( BadFast5Field &bf5 ){

					failed++;
					prog++;
					continue;
				}

				/*get the subsequence of the reference this read mapped to */
				r.referenceSeqMappedTo = reference.at(r.referenceMappedTo).substr(r.refStart, r.refEnd - r.refStart);

				/*pass on reads that cover a subsequences of the reference shorter than the minimum length */
				if ( r.referenceSeqMappedTo.size() < args.minL ) continue;

				//fetch the basecall from the bam file
				r.basecall = getQuerySequence(buffer[i]);

				//account for reverse complements
				if ( bam_is_rev(buffer[i]) ){

					r.basecall = reverseComplement( r.basecall );
					r.referenceSeqMappedTo = reverseComplement( r.referenceSeqMappedTo );
					r.isReverse = true;
				}

				normaliseEvents(r);

				//catch reads with rough event alignments that fail the QC
				if ( r.eventAlignment.size() == 0 ){

					failed++;
					prog++;
					continue;
				}

				std::stringstream ss; 
				llAcrossRead(r, windowLength, analogueModel, ss);

				#pragma omp critical
				{
					outFile << ss.rdbuf();
					prog++;
					pb.displayProgress( prog, failed );
				}
			}
			for ( unsigned int i = 0; i < buffer.size(); i++ ) bam_destroy1(buffer[i]);
			buffer.clear();
		}
		pb.displayProgress( prog, failed );	
	} while (result > 0);
	sam_itr_destroy(itr);
	std::cout << std::endl;
	return 0;
}
