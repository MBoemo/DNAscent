//----------------------------------------------------------
// Copyright 2019-2024 University of Oxford
// Written by Michael A. Boemo (mb915@cam.ac.uk)
// This software is licensed under GPL-3.0.  You should have
// received a copy of the license with this software.  If
// not, please Email the author.
//----------------------------------------------------------

#include <fstream>
#include <map>
#include <cmath>
#include <cstring>
#include <mutex>
#include <stdlib.h>
#include "refSignal.h"
#include "common.h"
#include "event_handling.h"
#include "alignment.h"
#include "pod5.h"
#include "fast5.h"
#include "htsInterface.h"
#include "data_IO.h"
#include "error_handling.h"
#include "config.h"
#include "../pod5-file-format/c++/pod5_format/c_api.h"
#include <omp.h>


static const char *help=
"refSignal: DNAscent executable that builds a reference nanopore signal profile across the genome.\n"
"To run DNAscent refSignal, do:\n"
"   DNAscent refSignal -b /path/to/sorted.bam -r /path/to/reference.fasta -i /path/to/index.dnascent -o /path/to/output\n"
"Required arguments are:\n"
"  -b,--bam                  path to sorted alignment BAM file,\n"
"  -r,--reference            path to genome reference in fasta format,\n"
"  -i,--index                path to DNAscent index,\n"
"  -o,--output               path to output file.\n"
"Optional arguments are:\n"
"  -t,--threads              number of threads (default is 1 thread),\n"
"  -q,--quality              minimum mapping quality (default is 20),\n"
"  -l,--length               minimum read length in bp (default is 1000).\n"
"DNAscent is under active development by the Boemo Group, Department of Pathology, University of Cambridge (https://www.boemogroup.org/).\n"
"Please submit bug reports to GitHub Issues (https://github.com/MBoemo/DNAscent/issues).";


struct RefSignalArguments {
	std::string bamFilename;
	std::string referenceFilename;
	std::string outputFilename;
	std::string indexFilename;
	int minQ;
	int minL;
	unsigned int threads;
};


RefSignalArguments parseRefSignalArguments( int argc, char** argv ){

	if ( argc < 2 ){

		std::cout << "Exiting with error.  Insufficient arguments passed to DNAscent refSignal." << std::endl << help << std::endl;
		exit(EXIT_FAILURE);
	}

	if ( std::string( argv[1] ) == "-h" or std::string( argv[1] ) == "--help" ){

		std::cout << help << std::endl;
		exit(EXIT_SUCCESS);
	}
	else if ( argc < 4 ){

		std::cout << "Exiting with error.  Insufficient arguments passed to DNAscent refSignal." << std::endl;
		exit(EXIT_FAILURE);
	}

	RefSignalArguments args;
	args.threads = 1;
	args.minQ    = 20;
	args.minL    = 1000;

	for ( int i = 1; i < argc; ){

		std::string flag( argv[i] );

		if ( flag == "-b" or flag == "--bam" ){

			if (i == argc-1) throw TrailingFlag(flag);
			args.bamFilename = std::string( argv[i+1] );
			i += 2;
		}
		else if ( flag == "-r" or flag == "--reference" ){

			if (i == argc-1) throw TrailingFlag(flag);
			args.referenceFilename = std::string( argv[i+1] );
			i += 2;
		}
		else if ( flag == "-t" or flag == "--threads" ){

			if (i == argc-1) throw TrailingFlag(flag);
			args.threads = std::stoi( argv[i+1] );
			i += 2;
		}
		else if ( flag == "-q" or flag == "--quality" ){

			if (i == argc-1) throw TrailingFlag(flag);
			args.minQ = std::stoi( argv[i+1] );
			if (args.minQ < 0) throw InvalidMappingThreshold();
			i += 2;
		}
		else if ( flag == "-l" or flag == "--length" ){

			if (i == argc-1) throw TrailingFlag(flag);
			args.minL = std::stoi( argv[i+1] );
			if (args.minL < 100) throw InvalidLengthThreshold();
			i += 2;
		}
		else if ( flag == "-i" or flag == "--index" ){

			if (i == argc-1) throw TrailingFlag(flag);
			args.indexFilename = std::string( argv[i+1] );
			i += 2;
		}
		else if ( flag == "-o" or flag == "--output" ){

			if (i == argc-1) throw TrailingFlag(flag);
			args.outputFilename = std::string( argv[i+1] );
			i += 2;
		}
		else throw InvalidOption( flag );
	}

	if (args.outputFilename == args.indexFilename or
	    args.outputFilename == args.referenceFilename or
	    args.outputFilename == args.bamFilename) throw OverwriteFailure();

	return args;
}


/*
 * Per-position accumulator using Welford's online algorithm so that
 * the full signal array never needs to be held in memory.
 * Each read contributes the mean of its aligned signal values at a
 * given reference position (one data point per read per position),
 * giving equal weight to every read regardless of dwell time.
 */
struct WelfordStats {
	long   count = 0;
	double mean  = 0.0;
	double M2    = 0.0;   // accumulated squared deviations from running mean

	void update( double x ){

		count++;
		double delta  = x - mean;
		mean         += delta / static_cast<double>(count);
		double delta2 = x - mean;
		M2           += delta * delta2;
	}

	double stddev() const {

		if (count < 2) return 0.0;
		return std::sqrt( M2 / static_cast<double>(count - 1) );
	}
};


int refSignal_main( int argc, char** argv ){

	RefSignalArguments args = parseRefSignalArguments( argc, argv );

	//load DNAscent index
	std::map< std::string, IndexEntry > readID2path;
	parseIndex( args.indexFilename, readID2path );

	//import fasta reference
	std::map< std::string, std::string > reference = import_reference_pfasta( args.referenceFilename );

	//open the bam
	std::cout << "Opening bam file... ";
	htsFile *bam_fh_cr = sam_open( args.bamFilename.c_str(), "r" );
	if (bam_fh_cr == NULL) throw IOerror( args.bamFilename );
	bam_hdr_t *bam_hdr_cr = sam_hdr_read( bam_fh_cr );
	std::cout << "ok." << std::endl;

	//open a log file
	std::cout << "Opening log file... ";
	std::string logFilename = strip_extension( args.outputFilename ) + ".refSignal.log";
	std::ofstream logfile( logFilename );
	if (logfile.is_open()) std::cout << "ok." << std::endl;
	else throw IOerror( logFilename );

	//initialise progress
	int numOfRecords = 0, prog = 0, failed = 0;
	countRecords( bam_fh_cr, bam_hdr_cr, numOfRecords, args.minQ, args.minL );
	progressBar pb( numOfRecords, true );
	bam_hdr_destroy( bam_hdr_cr );
	hts_close( bam_fh_cr );

	/*
	 * Genome-wide accumulator:  contig  ->  (reference position -> WelfordStats)
	 * Using a sparse map so only positions with actual coverage are stored.
	 */
	std::map< std::string, std::map< int, WelfordStats > > genomeStats;
	std::mutex stats_mtx;
	std::mutex log_mtx;

	pod5_init();

	int failedEvents = 0;
	unsigned int maxBufferSize;
	std::vector< bam1_t * > buffer;
	if ( args.threads <= 4 ) maxBufferSize = args.threads;
	else maxBufferSize = 4 * args.threads;

	htsFile   *bam_fh  = sam_open( args.bamFilename.c_str(), "r" );
	if (bam_fh == NULL) throw IOerror( args.bamFilename );
	bam_hdr_t *bam_hdr = sam_hdr_read( bam_fh );
	bam1_t    *itr_record = bam_init1();
	int result = sam_read1( bam_fh, bam_hdr, itr_record );

	while ( result >= 0 ){

		bam1_t *record = bam_dup1( itr_record );

		int mappingQual = record->core.qual;
		int refStart, refEnd;
		getRefEnd( record, refStart, refEnd );
		int queryLen = record->core.l_qseq;

		if ( mappingQual >= args.minQ and refEnd - refStart >= args.minL and queryLen != 0 ){
			buffer.push_back( record );
		}
		else{
			bam_destroy1( record );
		}

		result = sam_read1( bam_fh, bam_hdr, itr_record );

		if ( buffer.size() >= maxBufferSize or (buffer.size() > 0 and result == -1) ){

			#pragma omp parallel for schedule(dynamic) shared(buffer,Pore_Substrate_Config,args,prog,failed,genomeStats,stats_mtx,log_mtx) num_threads(args.threads)
			for ( unsigned int i = 0; i < buffer.size(); i++ ){

				DNAscent::read r( buffer[i], bam_hdr, readID2path, reference );

				if (r.refMismatch){

					std::cerr << std::endl << "Error: contig '" << r.referenceMappedTo
					          << "' found in BAM file but not in the reference genome." << std::endl;
					std::cerr << "Please check that the reference genome passed with -r matches the one used for alignment." << std::endl;
					exit(1);
				}

				if (r.missing){

					std::lock_guard<std::mutex> lock(log_mtx);
					logfile << "ReadID " << r.readID << " missing from index. Skipping." << std::endl;
					prog++;
					continue;
				}

				const char *ext = get_ext( r.filename.c_str() );

				if (strcmp(ext, "pod5") == 0){
					pod5_getSignal(r);
				}
				else if (strcmp(ext, "fast5") == 0){
					fast5_getSignal(r);
				}

				bool useFitPoreModel = false;
				normaliseEvents( r, useFitPoreModel );

				if ( r.eventAlignment.size() == 0 ){
					failed++;
					prog++;
					continue;
				}

				eventalign( r, Pore_Substrate_Config.windowLength_align );

				if ( not r.QCpassed ){
					failed++;
					prog++;
					continue;
				}

				/*
				 * Merge this read's per-position signal means into the genome-wide
				 * accumulator.  Each AlignedPosition may carry multiple raw signal
				 * samples from a single read; we reduce these to their mean so that
				 * every read contributes exactly one data point per covered position,
				 * regardless of dwell time.
				 */
				{
					std::lock_guard<std::mutex> lock(stats_mtx);
					std::map< int, WelfordStats > &contigStats = genomeStats[r.referenceMappedTo];
					for ( auto it = r.refCoordToAP.begin(); it != r.refCoordToAP.end(); ++it ){
						int pos = it->first;
						double sig_mean = it->second->getSignalMean();
						if ( sig_mean > 0.0 and sig_mean < 250.0 ){
							contigStats[pos].update( sig_mean );
						}
					}
				}

				prog++;
				pb.displayProgress( prog, failed, failedEvents );
			}
			buffer.clear();
		}
		pb.displayProgress( prog, failed, failedEvents );
	}

	bam_destroy1( itr_record );
	bam_hdr_destroy( bam_hdr );
	hts_close( bam_fh );
	std::cout << std::endl;
	pod5_terminate();

	//write the output table
	std::cout << "Writing output... ";
	std::ofstream outFile( args.outputFilename );
	if ( not outFile.is_open() ) throw IOerror( args.outputFilename );

	outFile << "#Alignment "      << args.bamFilename       << "\n";
	outFile << "#Genome "         << args.referenceFilename << "\n";
	outFile << "#Index "          << args.indexFilename     << "\n";
	outFile << "#Threads "        << args.threads           << "\n";
	outFile << "#MappingQuality " << args.minQ              << "\n";
	outFile << "#MappingLength "  << args.minL              << "\n";
	outFile << "#Version "        << VERSION                << "\n";
	outFile << "#Commit "         << getGitCommit()         << "\n";

	for ( auto cit = genomeStats.begin(); cit != genomeStats.end(); ++cit ){

		const std::string &contig = cit->first;
		const std::map< int, WelfordStats > &posMap = cit->second;

		for ( auto pit = posMap.begin(); pit != posMap.end(); ++pit ){

			const WelfordStats &s = pit->second;
			outFile << contig    << "\t"
			        << pit->first << "\t"
			        << s.mean    << "\t"
			        << s.stddev() << "\t"
			        << s.count   << "\n";
		}
	}

	outFile.close();
	std::cout << "ok." << std::endl;

	logfile.close();
	return 0;
}
