//----------------------------------------------------------
// Copyright 2019-2020 University of Oxford
// This software is licensed under GPL-3.0.  You should have
// received a copy of the license with this software.  If
// not, please Email the author.
//----------------------------------------------------------

//#define TEST_HMM 1
//#define TEST_LL 1
//#define TEST_ALIGNMENT 1
//#define TEST_METHYL 1

#include <fstream>
#include "detect.h"
#include <math.h>
#include <stdlib.h>
#include <limits>
#include "common.h"
#include "event_handling.h"
#include "probability.h"
#include "../fast5/include/fast5.hpp"
#include "../pod5-file-format/c++/pod5_format/c_api.h"
#include "detect.h"
#include "alignment.h"
#include "htsInterface.h"
#include "error_handling.h"
#include "tensor.h"
#include "config.h"
#include "pod5.h"
#include "fast5.h"
#include "reads.h"
#include "../htslib/htslib/hts.h"
#include "../htslib/htslib/sam.h"
#include <mutex>

static const char *help=
"trainCNN: DNAscent executable that generates HMM or CNN bootstrapped calls to build training data for DNAscent training.\n"
"Note: This executable is geared towards developers.\n"
"To run DNAscent trainCNN, do:\n"
"   DNAscent trainCNN -b /path/to/alignment.bam -r /path/to/reference.fasta -i /path/to/index.dnascent -o /path/to/output.trainCNN\n"
"Required arguments are:\n"
"  -b,--bam                  path to alignment BAM file,\n"
"  -r,--reference            path to genome reference in fasta format,\n"
"  -i,--index                path to DNAscent index,\n"
"  -o,--output               path to output file that will be generated.\n"
"Optional arguments are:\n"
"  -t,--threads              number of threads (default is 1 thread),\n"
"  --GPU                     use the GPU device indicated for prediction (default is CPU),\n"
"  -m,--maxReads             maximum number of reads to consider,\n"
"  -q,--quality              minimum mapping quality (default is 20),\n"
"  -l,--length               minimum read length in bp (default is 100),\n"
"     --HMM                  use HMM bootstrapping (default is CNN).\n"
"DNAscent is under active development by the Boemo Group, Department of Pathology, University of Cambridge (https://www.boemogroup.org/).\n"
"Please submit bug reports to GitHub Issues (https://github.com/MBoemo/DNAscent/issues).";

struct Arguments {
	std::string bamFilename;
	std::string referenceFilename;
	std::string outputFilename;
	std::string indexFilename;
	bool capReads;
	bool useHMM = false;
	bool useGPU = false;
	unsigned char GPUdevice = '0';
	int minQ, maxReads;
	int minL;
	unsigned int threads;
};

Arguments parseDataArguments( int argc, char** argv ){

	if( argc < 2 ){

		std::cout << "Exiting with error.  Insufficient arguments passed to DNAscent detect." << std::endl << help << std::endl;
		exit(EXIT_FAILURE);
	}

	if ( std::string( argv[ 1 ] ) == "-h" or std::string( argv[ 1 ] ) == "--help" ){

		std::cout << help << std::endl;
		exit(EXIT_SUCCESS);
	}
	else if( argc < 4 ){ //PLP&SY: check with Mike

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
			i+=2;
		}
		else if ( flag == "-l" or flag == "--length" ){

			if (i == argc-1) throw TrailingFlag(flag);

			std::string strArg( argv[ i + 1 ] );
			args.minL = std::stoi( strArg.c_str() );
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
			args.outputFilename = strArg;
			i+=2;
		}
		else if ( flag == "-m" or flag == "--maxReads" ){

			if (i == argc-1) throw TrailingFlag(flag);

			std::string strArg( argv[ i + 1 ] );
			args.capReads = true;
			args.maxReads = std::stoi( strArg.c_str() );
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



int data_main( int argc, char** argv ){

	Arguments args = parseDataArguments( argc, argv );

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

	std::ofstream outFile( args.outputFilename );
	if ( not outFile.is_open() ) throw IOerror( args.outputFilename );

	//load the bam
	std::cout << "Opening bam file... ";
	htsFile *bam_fh_cr = sam_open((args.bamFilename).c_str(), "r");
	if (bam_fh_cr == NULL) throw IOerror(args.bamFilename);
	bam_hdr_t *bam_hdr_cr = sam_hdr_read(bam_fh_cr);
	std::cout << "ok." << std::endl;

	//open a log file
	std::cout << "Opening log file... ";
	std::string logFilename = strip_extension(args.outputFilename);
	logFilename += ".trainCNN.log";
	std::ofstream logfile(logFilename);
	if (logfile.is_open()) std::cout << "ok." << std::endl;
	else throw IOerror(logFilename);
	std::mutex mtx;

	//initialise progress
	int numOfRecords = 0, prog = 0, failed = 0;
	countRecords( bam_fh_cr, bam_hdr_cr, numOfRecords, args.minQ, args.minL );
	progressBar pb(numOfRecords,true);

	pod5_init();

	int failedEvents = 0;
	unsigned int maxBufferSize;
	std::vector< bam1_t * > buffer;
	if ( args.threads <= 4 ) maxBufferSize = args.threads; //PLP&SY: check with Mike
	else maxBufferSize = 4*(args.threads);

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

			#pragma omp parallel for schedule(dynamic) shared(buffer,Pore_Substrate_Config,args,prog,failed,session,inputOps) num_threads(args.threads)
			for (unsigned int i = 0; i < buffer.size(); i++){

				DNAscent::read r(buffer[i], bam_hdr, readID2path, reference);

				if (r.missing){

					std::cerr << "ReadID " << r.readID << " missing from index. Skipping." << std::endl;
					prog++;
					continue;
				}

				const char *ext = get_ext(r.filename.c_str());
				
				if (strcmp(ext,"pod5") == 0){
					pod5_getSignal(r);
				}
				else if (strcmp(ext,"fast5") == 0){
					fast5_getSignal(r);
				}

				bool useFitPoreModel = false;	
				normaliseEvents(r, useFitPoreModel);							
				if ( (r.eventAlignment).size() == 0 ){

					failed++;
					prog++;
					continue;
				}
				
				//do a first event alignment to make DNN input tensors
				eventalign( r, Pore_Substrate_Config.windowLength_align);
				
				//run analogue prediction
				if (args.useHMM) llAcrossRead(r, 12);
				else runCNN(r,session,inputOps, true);
				
				//clear the read and re-annotate wtih the DNN analogue predictions
				eventalign(r, Pore_Substrate_Config.windowLength_align);				

				if (not r.QCpassed){
					failed++;
					prog++;
					continue;
				}

				#pragma omp critical
				{
					outFile << r.humanReadable_eventalignOut;
					prog++;
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
	std::cout << std::endl;
	pod5_terminate();
	return 0;
}
