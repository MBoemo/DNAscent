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
#include "alignment.h"
#include "pod5.h"
#include "fast5.h"
#include "probability.h"
#include "../fast5/include/fast5.hpp"
#include "../tensorflow/include/tensorflow/c/eager/c_api.h"
#include "../pod5-file-format/c++/pod5_format/c_api.h"
#include "htsInterface.h"
#include "error_handling.h"
#include "config.h"
#include <omp.h>


static const char *help=
"detect: DNAscent executable that detects BrdU and EdU in Oxford Nanopore reads.\n"
"To run DNAscent detect, do:\n"
"   DNAscent detect -b /path/to/alignment.bam -r /path/to/reference.fasta -i /path/to/index.dnascent -o /path/to/output.detect\n"
"Required arguments are:\n"
"  -b,--bam                  path to alignment BAM file,\n"
"  -r,--reference            path to genome reference in fasta format,\n"
"  -i,--index                path to DNAscent index,\n"
"  -o,--output               path to output file with extension `detect` (for human-readable format) or `bam` (for modbam format).\n"
"Optional arguments are:\n"
"  -t,--threads              number of threads (default is 1 thread),\n"
"  --GPU                     use the GPU device indicated for prediction (default is CPU),\n"
"  -q,--quality              minimum mapping quality (default is 20),\n"
"  -l,--length               minimum read length in bp (default is 1000).\n"
"DNAscent is under active development by the Boemo Group, Department of Pathology, University of Cambridge (https://www.boemogroup.org/).\n"
"Please submit bug reports to GitHub Issues (https://github.com/MBoemo/DNAscent/issues).";

struct Arguments {
	std::string bamFilename;
	std::string referenceFilename;
	std::string outputFilename;
	std::string indexFilename;
	bool humanReadable = false;
	bool useGPU = false;
	bool useHMM = false;
	unsigned char GPUdevice = '0';
	int minQ = 20;
	int minL = 1000;
	unsigned int threads = 1;
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

			if (i == argc-1) throw TrailingFlag(flag);

			std::string strArg( argv[ i + 1 ] );
			args.bamFilename = strArg;
			i+=2;
		}
		else if ( flag == "-r" or flag == "--reference" ){

			if (i == argc-1) throw TrailingFlag(flag);

			std::string strArg( argv[ i + 1 ] );
			args.referenceFilename = strArg;
			i+=2;
		}
		else if ( flag == "-t" or flag == "--threads" ){

			if (i == argc-1) throw TrailingFlag(flag);

			std::string strArg( argv[ i + 1 ] );
			args.threads = std::stoi( strArg.c_str() );
			i+=2;
		}
		else if ( flag == "-q" or flag == "--quality" ){

			if (i == argc-1) throw TrailingFlag(flag);

			std::string strArg( argv[ i + 1 ] );
			args.minQ = std::stoi( strArg.c_str() );
			
			if (args.minQ < 0) throw InvalidMappingThreshold();
			
			i+=2;
		}
		else if ( flag == "-l" or flag == "--length" ){

			if (i == argc-1) throw TrailingFlag(flag);

			std::string strArg( argv[ i + 1 ] );
			args.minL = std::stoi( strArg.c_str() );
			
			if (args.minL < 100) throw InvalidLengthThreshold();
			if (args.minL < 1000) std::cerr << "Warning: DNAscent detect may show inaccuracies or high fail rates on short reads (< 1 kb)." << std::endl;
			
			i+=2;
		}
		else if ( flag == "-i" or flag == "--index" ){

			if (i == argc-1) throw TrailingFlag(flag);

			std::string strArg( argv[ i + 1 ] );
			args.indexFilename = strArg;
			i+=2;
		}
		else if ( flag == "-o" or flag == "--output" ){

			if (i == argc-1) throw TrailingFlag(flag);

			std::string strArg( argv[ i + 1 ] );
			
			const char *ext = get_ext(strArg.c_str());
				
			if (strcmp(ext,"bam") == 0){
				args.humanReadable = false;
			}
			else if (strcmp(ext,"detect") == 0){
				args.humanReadable = true;
			}
			else{
				throw InvalidExtension(ext);
			}
			
			args.outputFilename = strArg;
			i+=2;
		}
		else if ( flag == "--HMM" ){
		
			args.useHMM = true;
			i+=1;
		}
		else if ( flag == "--GPU" ){
		
			if (i == argc-1) throw TrailingFlag(flag);		
		
			args.useGPU = true;
			std::string strArg( argv[ i + 1 ] );
			if (strArg.length() > 1) throw InvalidDevice(strArg);

			args.GPUdevice = *argv[ i + 1 ];

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

	unsigned int k = Pore_Substrate_Config.kmer_len;

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

	std::vector< double > I_curr(2*windowSize, NAN), D_curr(2*windowSize, NAN), M_curr(2*windowSize, NAN), I_prev(2*windowSize, NAN), D_prev(2*windowSize, NAN), M_prev(2*windowSize, NAN);
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

		std::string kmer = sequence.substr(0, k);

		std::pair<double,double> meanStd = Pore_Substrate_Config.unlabelled_model[kmer2index(kmer, k)];
		
		level_mu = meanStd.first;
		level_sigma = meanStd.second;

		matchProb = eln( normalPDF( level_mu, level_sigma, (observations[t] - scalings.shift)/scalings.scale ) );
		insProb = 0.0; //log(1) = 0; set probability equal to 1 and use transition probability as weighting

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
			kmer = sequence.substr(i, k);
			insProb = 0.0; //log(1) = 0; set probability equal to 1 and use transition probability as weighting
			if ( useBrdU and BrdUStart <= i and i <= BrdUEnd and kmer.find('T') != std::string::npos ){
				
				std::pair<double,double> analogue_meanStd = Pore_Substrate_Config.analogue_model[kmer2index(kmer, k)];
				level_mu = analogue_meanStd.first;
				level_sigma = analogue_meanStd.second;
				matchProb = eln( normalPDF( level_mu, level_sigma, (observations[t] - scalings.shift)/scalings.scale ) );
			}
			else{

				std::pair<double,double> meanStd = Pore_Substrate_Config.unlabelled_model[kmer2index(kmer, k)];
				level_mu = meanStd.first;
				level_sigma = meanStd.second;
				matchProb = eln( normalPDF( level_mu, level_sigma, (observations[t] - scalings.shift)/scalings.scale ) );
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
std::cerr << scalings.shift << " " << scalings.scale << std::endl;
std::cerr << sequence << std::endl;
for (auto ob = observations.begin(); ob < observations.end(); ob++){
	std::cerr << *ob << " ";
}
std::cerr << std::endl;
std::cerr << forwardProb << std::endl;
#endif

	return forwardProb;
}


std::vector< unsigned int > getPOIs( std::string &refSeq, int windowLength ){

	std::vector< unsigned int > POIs;

	for ( unsigned int i = 2*windowLength; i < refSeq.length() - 2*windowLength; i++ ){

		if (refSeq.substr(i,1) == "T") POIs.push_back(i);
	}
	return POIs;
}


HMMdetection llAcrossRead( DNAscent::read &r, unsigned int windowLength){
                          
	HMMdetection hmm_out;

	std::map<unsigned int, double> refPos2likelihood;

	unsigned int k = Pore_Substrate_Config.kmer_len;

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

	hmm_out.stdout += ">" + r.readID + " " + r.referenceMappedTo + " " + std::to_string(r.refStart) + " " + std::to_string(r.refEnd) + " " + strand + "\n";

	for ( unsigned int i = 0; i < POIs.size(); i++ ){

		unsigned int posOnRef = POIs[i];
		unsigned int posOnQuery = (r.refToQuery).at(posOnRef);

		std::string readSnippet = (r.referenceSeqMappedTo).substr(posOnRef - windowLength, 2*windowLength + k);

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
		unsigned int spanOnQuery = (r.refToQuery)[posOnRef + windowLength] - (r.refToQuery)[posOnRef - windowLength];
		assert(spanOnQuery >= 0);

		/*get the events that correspond to the read snippet */
		bool first = true;
		if ( r.isReverse ){

			for ( int j = readHead; j >= 0; j-- ){

				assert(j >= 0);
				assert(posOnRef - windowLength >= 0);
				assert(posOnRef + windowLength < (r.refToQuery).size());

				/*if an event has been aligned to a position in the window, add it */
				if ( (r.refToQuery)[posOnRef - windowLength] <= (r.eventAlignment)[j].second and (r.eventAlignment)[j].second < (r.refToQuery)[posOnRef + windowLength] ){

					if (first){
						readHead = j;
						first = false;
						//std::cout << "READHEAD:" << j << " " << readHead << std::endl;
					}
					
					double ev = (r.events)[(r.eventAlignment)[j].first].mean;
					if (ev > 0. and ev < 250.0){

						eventSnippet.push_back(ev);
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

				assert(posOnRef - windowLength >= 0);
				assert(posOnRef + windowLength < (r.refToQuery).size());				

				/*if an event has been aligned to a position in the window, add it */
				if ( (r.refToQuery)[posOnRef - windowLength] <= (r.eventAlignment)[j].second and (r.eventAlignment)[j].second < (r.refToQuery)[posOnRef + windowLength] ){

					if (first){
						readHead = j;
						first = false;
						//std::cout << "READHEAD:" << j << " " << readHead << std::endl;
					}
					
					double ev = (r.events)[(r.eventAlignment)[j].first].mean;
					if (ev > 0. and ev < 250.0){
						eventSnippet.push_back(ev);
					}
				}

				/*stop once we get to the end of the window */
				if ( (r.eventAlignment)[j].second >= (r.refToQuery)[posOnRef + windowLength] ) break;
			}
		}

		//catch abnormally few events
		if ( eventSnippet.size() < 2*windowLength - k ) continue;

		/*
		TESTING - print out the read snippet, the ONT model, and the aligned events
		std::cout << readSnippet << std::endl;
		for ( int pos = 0; pos < readSnippet.length()-5; pos++ ){

			std::cout << readSnippet.substr(pos,6) << "\t" << Pore_Substrate_Config.pore_model.at( readSnippet.substr(pos,6) ).first << std::endl;
		}
		for ( auto ev = eventSnippet.begin(); ev < eventSnippet.end(); ev++){
			double scaledEv =  (*ev - r.scalings.shift) / r.scalings.scale;
			std::cout << scaledEv << std::endl;
		}
		*/

		//calculate where we are on the assembly - if we're a reverse complement, we're moving backwards down the reference genome
		int globalPosOnRef;
		assert(posOnQuery < (r.basecall).size());
		std::string kmerQuery = (r.basecall).substr(posOnQuery, k);
		assert(posOnRef < (r.referenceSeqMappedTo).size());		
		std::string kmerRef = (r.referenceSeqMappedTo).substr(posOnRef, k);
		if ( r.isReverse ){

			globalPosOnRef = r.refEnd - posOnRef - 1;
			kmerQuery = reverseComplement( kmerQuery );
			kmerRef = reverseComplement( kmerRef );
		}
		else{

			globalPosOnRef = r.refStart + posOnRef;
		}

		//make the BrdU call
		std::string kOI = (r.referenceSeqMappedTo).substr(posOnRef-4,k);
		size_t BrdUStart = windowLength - (int) k/2;
		size_t BrdUEnd = windowLength + (int) k/2;
		double logProbAnalogue = sequenceProbability( eventSnippet, readSnippet, windowLength, true, r.scalings, BrdUStart, BrdUEnd );
		double logProbThymidine = sequenceProbability( eventSnippet, readSnippet, windowLength, false, r.scalings, 0, 0 );
		double logLikelihoodRatio = logProbAnalogue - logProbThymidine;


#if TEST_LL
double runningKL = 0.0;
for (unsigned int s = 0; s < readSnippet.length() - k; s++){
	std::string kmer = readSnippet.substr(s,Pore_Substrate_Config.kmer_len);
	if ( BrdUStart <= s and s <= BrdUEnd and kmer.find('T') != std::string::npos and Pore_Substrate_Config.analogue_model.count(kmer) > 0 ){
		runningKL += KLdivergence( Pore_Substrate_Config.pore_model.at(kmer).first, Pore_Substrate_Config.pore_model.at(kmer).second, Pore_Substrate_Config.analogue_model.at(kmer).first, Pore_Substrate_Config.analogue_model.at(kmer).second );
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

		hmm_out.stdout += std::to_string(globalPosOnRef) + "\t" + std::to_string(logLikelihoodRatio) + "\t" + kmerRef + "\t" + kmerQuery + "\n";

		//adjust reference position so that we make the call at the middle of the kmer
		if ( r.isReverse ){

			globalPosOnRef += (int) k/2;
		}
		else{

			globalPosOnRef -= (int) k/2;
		}
		
		hmm_out.refposToLikelihood[globalPosOnRef] = std::make_pair(logLikelihoodRatio,0.);
	}
	return hmm_out;
}


void runCNN(DNAscent::read &r, std::shared_ptr<ModelSession> session, std::vector<TF_Output> inputOps, bool humanReadable){

	int NumInputs = 3;
	int NumOutputs = 1;

	std::vector<TF_Tensor*> input_tensors;
	TF_Tensor* OutputValues;

	//core sequence input
	std::vector<size_t> protoSequenceShape = r.getSequenceShape();
	TensorShape input_sequenceShape={{1, (int64_t) protoSequenceShape[0]}, 2};

	std::vector<float> unformattedCoreSequenceTensor = r.makeCoreSequenceTensor();

	size_t sizeSequence = unformattedCoreSequenceTensor.size();
	assert(sizeSequence > 0);

	float *tmp_coreSequenceArray = (float *)malloc(sizeSequence*sizeof(float));
	for(size_t i = 0; i < sizeSequence; i++){
		tmp_coreSequenceArray[i] = unformattedCoreSequenceTensor[i];
	}

	TF_Tensor* CoreSequenceInputTensor = TF_NewTensor(TF_FLOAT,
		input_sequenceShape.values,
		input_sequenceShape.dim,
		(void *)tmp_coreSequenceArray,
		sizeSequence*sizeof(float),
		cpp_array_deallocator<float>,
		nullptr);

	input_tensors.push_back(CoreSequenceInputTensor);
	
	//residual sequence input (inherits the same shape as core sequence)
	std::vector<float> unformattedResidualSequenceTensor = r.makeResidualSequenceTensor();

	float *tmp_resSequenceArray = (float *)malloc(sizeSequence*sizeof(float));
	for(size_t i = 0; i < sizeSequence; i++){
		tmp_resSequenceArray[i] = unformattedResidualSequenceTensor[i];
	}

	TF_Tensor* ResidualSequenceInputTensor = TF_NewTensor(TF_FLOAT,
		input_sequenceShape.values,
		input_sequenceShape.dim,
		(void *)tmp_resSequenceArray,
		sizeSequence*sizeof(float),
		cpp_array_deallocator<float>,
		nullptr);

	input_tensors.push_back(ResidualSequenceInputTensor);

	//signal input
	std::vector<size_t> protoSignalShape = r.getSignalShape();
	TensorShape input_signalShape={{1, (int64_t) protoSignalShape[0], (int64_t) protoSignalShape[1], (int64_t) protoSignalShape[2]}, 4};

	std::vector<float> unformattedSignalTensor = r.makeSignalTensor();

	size_t sizeSignal = unformattedSignalTensor.size();
	assert(sizeSignal > 0);

	float *tmp_signalArray = (float *)malloc(sizeSignal*sizeof(float));
	for(size_t i = 0; i < sizeSignal; i++){
		tmp_signalArray[i] = unformattedSignalTensor[i];
	}

	TF_Tensor* SignalInputTensor = TF_NewTensor(TF_FLOAT,
		input_signalShape.values,
		input_signalShape.dim,
		(void *)tmp_signalArray,
		sizeSignal*sizeof(float),
		cpp_array_deallocator<float>,
		nullptr);

	input_tensors.push_back(SignalInputTensor);

	//Run the Session
	CStatus status;
	TF_SessionRun(*(session->session.get()), NULL, &inputOps[0], &input_tensors[0], NumInputs, &session->outputs, &OutputValues, NumOutputs, NULL, 0, NULL, status.ptr);

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
	std::vector<unsigned int> refCoordinates = r.getReferenceCoords();
	std::vector<unsigned int> refIndices = r.getReferenceIndices();
	std::vector<unsigned int> queryIndices = r.getQueryIndices();
	std::vector<std::string> kmers = r.getKmers();

	size_t output_size = TF_TensorByteSize(OutputValues) / sizeof(float);
	assert(output_size == protoSequenceShape[0] * outputFields);
	float *output_array = (float *)TF_TensorData(OutputValues);

	//write the output
	unsigned int pos_ctr = 0;
	std::vector<std::string> lines;
	lines.reserve(refCoordinates.size());
	unsigned int thisRefCoord = refCoordinates[0];
	std::string str_line;
	
	if (humanReadable) r.humanReadable_detectOut += ">" + r.readID + " " + r.referenceMappedTo + " " + std::to_string(r.refStart) + " " + std::to_string(r.refEnd) + " " + r.strand + "\n";
	
	for(size_t i = 0; i < output_size; i++){
		if((i+1)%outputFields==0){

			//only output kmers where the middle position is a T
			if (kmers[pos_ctr].substr(4,1) != "T"){ 
				pos_ctr++;
				continue;
			}
			
			r.refCoordToCalls[thisRefCoord] = std::make_pair(output_array[i], output_array[i-1]);

			if (humanReadable){
				str_line += std::to_string(thisRefCoord) + "\t" + std::to_string(output_array[i])+ "\t" + std::to_string(output_array[i-1]);
				if (r.isReverse) str_line += "\t" + reverseComplement(kmers[pos_ctr]);
				else str_line += "\t" + kmers[pos_ctr];
				lines.push_back(str_line);
				str_line = "";
			}
			else if ( not r.refToDel.at(refIndices[pos_ctr]) ){
			
				r.queryIndexToCalls[queryIndices[pos_ctr]] = std::make_pair(output_array[i], output_array[i-1]);			
			}
			
			pos_ctr++;
		}
		else{
			if (i != output_size-1) thisRefCoord = refCoordinates[pos_ctr];
		}
	}

	TF_DeleteTensor(OutputValues);
	TF_DeleteTensor(CoreSequenceInputTensor);
	TF_DeleteTensor(ResidualSequenceInputTensor);
	TF_DeleteTensor(SignalInputTensor);

	if (humanReadable){
		if (r.isReverse) std::reverse(lines.begin(),lines.end());

		for (auto s = lines.begin(); s < lines.end(); s++){
			r.humanReadable_detectOut += *s + "\n";
		}
	}
	else{
	
		r.writeModBamTag();
	}
}


int detect_main( int argc, char** argv ){

	Arguments args = parseDetectArguments( argc, argv );

	//load DNAscent index
	std::map< std::string, IndexEntry > readID2path;
	parseIndex( args.indexFilename, readID2path );

	//get the neural network model path
	std::string pathExe = getExePath();
	std::string modelPath = pathExe + Pore_Substrate_Config.fn_dnn_model;
	std::string input1_layer_name = Pore_Substrate_Config.dnn_model_inputLayer1;
	std::string input2_layer_name = Pore_Substrate_Config.dnn_model_inputLayer2;
	std::string input3_layer_name = Pore_Substrate_Config.dnn_model_inputLayer3;

	std::pair< std::shared_ptr<ModelSession>, std::shared_ptr<TF_Graph *> > modelPair;

	if (not args.useGPU){

		modelPair = model_load_cpu_twoInputs(modelPath.c_str(), args.threads);
	}
	else{

		modelPair = model_load_gpu_twoInputs(modelPath.c_str(), args.GPUdevice, args.threads);
	}

	std::shared_ptr<ModelSession> session = modelPair.first;
	std::shared_ptr<TF_Graph *> Graph = modelPair.second;

	auto input1_op = TF_GraphOperationByName(*(Graph.get()), input1_layer_name.c_str());
	auto input2_op = TF_GraphOperationByName(*(Graph.get()), input2_layer_name.c_str());
	auto input3_op = TF_GraphOperationByName(*(Graph.get()), input3_layer_name.c_str());
	if(!input1_op or !input2_op or !input3_op){
		std::cout << "bad input name" << std::endl;
		exit(0);
	}

	std::vector<TF_Output> inputOps = {{input1_op,0}, {input2_op,0}, {input3_op,0}};

	//import fasta reference
	std::map< std::string, std::string > reference = import_reference_pfasta( args.referenceFilename );

	//load the bam
	std::cout << "Opening bam file... ";
	htsFile *bam_fh_cr = sam_open((args.bamFilename).c_str(), "r");
	if (bam_fh_cr == NULL) throw IOerror(args.bamFilename);
	bam_hdr_t *bam_hdr_cr = sam_hdr_read(bam_fh_cr);
	std::cout << "ok." << std::endl;

	//make the output writer
	OutputFormat format;
	if (args.humanReadable) format = OutputFormat::HumanReadable;
	else format = OutputFormat::Sam;
	auto writer = OutputWriterFactory::createWriter(format);
	writer -> open(args.outputFilename);

	//write the output header
	if (args.humanReadable){
		std::string outHeader = writeDetectHeader(args.bamFilename, args.referenceFilename, args.indexFilename, args.threads, false, args.minQ, args.minL, args.useGPU);	
		writer -> writeHeader_HR(outHeader);
	}
	else{
		writer -> writeHeader_sam(bam_hdr_cr);
	}

	//initialise progress
	int numOfRecords = 0, prog = 0, failed = 0;
	countRecords( bam_fh_cr, bam_hdr_cr, numOfRecords, args.minQ, args.minL );
	progressBar pb(numOfRecords,true);
	bam_hdr_destroy(bam_hdr_cr);
	hts_close(bam_fh_cr);

	pod5_init();

	int failedEvents = 0;
	unsigned int maxBufferSize;
	std::vector< bam1_t * > buffer;
	maxBufferSize = 16*(args.threads);

	htsFile *bam_fh = sam_open((args.bamFilename).c_str(), "r");
	if (bam_fh == NULL) throw IOerror(args.bamFilename);
	bam_hdr_t *bam_hdr = sam_hdr_read(bam_fh);
	bam1_t *itr_record = bam_init1();
	int result = sam_read1(bam_fh, bam_hdr, itr_record);

	while(result >= 0){
	
		bam1_t *record = bam_dup1(itr_record);

		//add the record to the buffer if it passes the user's criteria, otherwise destroy it cleanly
		int mappingQual = record -> core.qual;
		int refStart,refEnd;
		getRefEnd(record,refStart,refEnd);
		int queryLen = record -> core.l_qseq;

		if ( mappingQual >= args.minQ and refEnd - refStart >= args.minL and queryLen != 0 ){

			buffer.push_back(record);
		}
		else{
			bam_destroy1(record);
		}

		result = sam_read1(bam_fh, bam_hdr, itr_record);

		//if we've filled up the buffer with reads, compute them in parallel
		if (buffer.size() >= maxBufferSize or (buffer.size() > 0 and result == -1 ) ){

			#pragma omp parallel for schedule(dynamic) shared(buffer,Pore_Substrate_Config,args,prog,failed,session,inputOps,writer) num_threads(args.threads)
			for (unsigned int i = 0; i < buffer.size(); i++){

				DNAscent::read r(buffer[i], bam_hdr, readID2path, reference);

				const char *ext = get_ext(r.filename.c_str());

				if (strcmp(ext,"pod5") == 0){
					pod5_getSignal(r);
				}
				else if (strcmp(ext,"fast5") == 0){
					fast5_getSignal(r);
				}

				//for HMM
				//bool useFitPoreModel = true;
				//normaliseEvents(r, useFitPoreModel);

				//for DNN
				bool useFitPoreModel = false;
				normaliseEvents( r, useFitPoreModel);

				//catch reads with rough event alignments that fail the QC
				if ( r.eventAlignment.size() == 0 ){
					failed++;
					prog++;
					continue;
				}

				//HMMdetection hmm_likelihood = llAcrossRead(r, 12);
				//readOut = hmm_likelihood.stdout;
				
				eventalign( r, Pore_Substrate_Config.windowLength_align);

				if (not r.QCpassed){
					failed++;
					prog++;
					continue;
				}

				runCNN(r,session,inputOps,args.humanReadable);

				prog++;
				pb.displayProgress( prog, failed, failedEvents );
				
				#pragma omp critical
				{
					writer -> write(r);
					pb.displayProgress( prog, failed, failedEvents );
				}
			}
			buffer.clear();
		}
		pb.displayProgress( prog, failed, failedEvents );
	}
	bam_destroy1(itr_record);
	bam_hdr_destroy(bam_hdr);
	hts_close(bam_fh);
	writer -> close();
	std::cout << std::endl;
	pod5_terminate();
	return 0;
}
