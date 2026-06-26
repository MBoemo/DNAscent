//----------------------------------------------------------
// Copyright 2019-2024 University of Oxford
// Written by Michael A. Boemo (mb915@cam.ac.uk)
// This software is licensed under GPL-3.0.  You should have
// received a copy of the license with this software.  If
// not, please Email the author.
//----------------------------------------------------------

#include <fstream>
#include <sstream>
#include <map>
#include <math.h>
#include <cstring>
#include <mutex>
#include <stdlib.h>
#include "refScan.h"
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
"refScan: DNAscent executable that scores each read against a reference signal\n"
"profile to identify candidate modification sites.\n"
"\n"
"For each read a rolling k-mer window (k = kmer length of the pore model, e.g. 9\n"
"for R10.4.1 or 6 for R9) is slid along the event alignment.  At each window\n"
"centre the observed per-position signal means are compared to the corresponding\n"
"Gaussian distributions stored in the reference signal file, and the mean\n"
"per-position log-likelihood is reported.  Positions where the score is more\n"
"negative than expected indicate signal deviations consistent with a modification.\n"
"\n"
"To run DNAscent refScan, do:\n"
"   DNAscent refScan -b /path/to/alignment.bam -r /path/to/reference.fasta \\\n"
"                         -i /path/to/index.dnascent -s /path/to/refSignal \\\n"
"                         -o /path/to/output\n"
"Required arguments are:\n"
"  -b,--bam                  path to alignment BAM file,\n"
"  -r,--reference            path to genome reference in fasta format,\n"
"  -i,--index                path to DNAscent index,\n"
"  -s,--refSignal            path to reference signal file produced by DNAscent refSignal,\n"
"  -o,--output               path to output file.\n"
"Optional arguments are:\n"
"  -t,--threads              number of threads (default is 1 thread),\n"
"  -q,--quality              minimum mapping quality (default is 20),\n"
"  -l,--length               minimum read length in bp (default is 1000),\n"
"  -w,--minWindowCov         minimum number of positions in the k-mer window that\n"
"                            must have both read signal and reference signal for a\n"
"                            window to be scored (default is 5).\n"
"                            Maximum is the k-mer length (9 for R10.4.1, 6 for R9).\n"
"DNAscent is under active development by the Boemo Group, Department of Pathology, University of Cambridge (https://www.boemogroup.org/).\n"
"Please submit bug reports to GitHub Issues (https://github.com/MBoemo/DNAscent/issues).";


struct RefScanArguments {
	std::string bamFilename;
	std::string referenceFilename;
	std::string refSignalFilename;
	std::string outputFilename;
	std::string indexFilename;
	int minQ;
	int minL;
	int minWindowCov;
	unsigned int threads;
};


RefScanArguments parseRefScanArguments( int argc, char** argv ){

	if ( argc < 2 ){

		std::cout << "Exiting with error.  Insufficient arguments passed to DNAscent refScan." << std::endl << help << std::endl;
		exit(EXIT_FAILURE);
	}

	if ( std::string( argv[1] ) == "-h" or std::string( argv[1] ) == "--help" ){

		std::cout << help << std::endl;
		exit(EXIT_SUCCESS);
	}
	else if ( argc < 4 ){

		std::cout << "Exiting with error.  Insufficient arguments passed to DNAscent refScan." << std::endl;
		exit(EXIT_FAILURE);
	}

	RefScanArguments args;
	args.threads      = 1;
	args.minQ         = 20;
	args.minL         = 1000;
	args.minWindowCov = 5;

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
		else if ( flag == "-s" or flag == "--refSignal" ){

			if (i == argc-1) throw TrailingFlag(flag);
			args.refSignalFilename = std::string( argv[i+1] );
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
		else if ( flag == "-w" or flag == "--minWindowCov" ){

			if (i == argc-1) throw TrailingFlag(flag);
			args.minWindowCov = std::stoi( argv[i+1] );
			int klen = static_cast<int>( Pore_Substrate_Config.kmer_len );
			if (args.minWindowCov < 1 or args.minWindowCov > klen){
				std::cout << "Exiting with error.  -w must be between 1 and " << klen << "." << std::endl;
				exit(EXIT_FAILURE);
			}
			i += 2;
		}
		else throw InvalidOption( flag );
	}

	if (args.outputFilename == args.indexFilename or
	    args.outputFilename == args.referenceFilename or
	    args.outputFilename == args.bamFilename or
	    args.outputFilename == args.refSignalFilename) throw OverwriteFailure();

	return args;
}


/*
 * One entry per reference position loaded from the refSignal output file.
 */
struct RefSigEntry {
	double mean;
	double std;
	long   coverage;
};


/*
 * Load the reference signal file produced by DNAscent refSignal.
 * Format:
 *   comment lines start with '#'
 *   first non-comment line is the header
 *   subsequent lines: contig \t position \t mean_signal \t std_signal \t coverage
 */
std::map< std::string, std::map< int, RefSigEntry > > loadRefSignal( const std::string &filename ){

	std::ifstream f( filename );
	if ( not f.is_open() ) throw IOerror( filename );

	std::map< std::string, std::map< int, RefSigEntry > > data;

	std::string line;
	while ( std::getline(f, line) ){

		if ( line.empty() ) continue;
		if ( line[0] == '#' ) continue;

		std::istringstream ss(line);
		std::string contig;
		int pos;
		RefSigEntry e;
		ss >> contig >> pos >> e.mean >> e.std >> e.coverage;
		data[contig][pos] = e;
	}

	return data;
}


/*
 * Score a read against the reference signal using a rolling k-mer window
 * (k = Pore_Substrate_Config.kmer_len; 9 for R10.4.1, 6 for R9).
 *
 * For each position in the read's event alignment that has signal, we treat
 * it as the centre of a k-position window.  For every position in
 * the window that has both an observed signal mean and a reference Gaussian
 * (mean + std from refSignal), we accumulate the per-position log-likelihood:
 *
 *   log N(obs; mu, sigma) = -0.5 * [log(2*pi) + 2*log(sigma) + ((obs-mu)/sigma)^2]
 *
 * The window score is the mean log-likelihood over valid positions.
 * Lower (more negative) scores indicate deviation from the reference and are
 * consistent with a modification at or near the window centre.
 *
 * The output string follows the same per-read header convention used by
 * DNAscent align:
 *   >readID  contig  refStart  refEnd  strand
 *   pos  meanLL  windowCoverage
 */
void scoreWindowsOnRead( DNAscent::read &r,
                         const std::map< std::string, std::map< int, RefSigEntry > > &refSig,
                         int minWindowCov,
                         std::string &out ){

	static const double LOG2PI = log( 2.0 * M_PI );

	// Window spans exactly kmer_len positions centred on each covered base.
	// For odd k (e.g. 9):  [center - k/2 .. center + k/2]      (9 positions)
	// For even k (e.g. 6): [center - k/2 .. center + k/2 - 1]  (6 positions)
	const int k      = static_cast<int>( Pore_Substrate_Config.kmer_len );
	const int half_l = k / 2;           // left extent  (inclusive)
	const int half_r = k - 1 - half_l; // right extent (inclusive)

	// Nothing to do if there is no reference signal for this contig
	auto cit = refSig.find( r.referenceMappedTo );
	if ( cit == refSig.end() ) return;
	const std::map< int, RefSigEntry > &contigRef = cit->second;

	out += ">" + r.readID
	     + " " + r.referenceMappedTo
	     + " " + std::to_string( r.refStart )
	     + " " + std::to_string( r.refEnd )
	     + " " + r.strand + "\n";

	// Iterate over every position for which this read has aligned signal and
	// use it as the centre of a k-mer window
	for ( auto ctr_it = r.refCoordToAP.begin(); ctr_it != r.refCoordToAP.end(); ++ctr_it ){

		int center = ctr_it->first;

		double window_ll = 0.0;
		int    valid     = 0;

		for ( int off = -half_l; off <= half_r; off++ ){

			int pos = center + off;

			// Observed signal at this position in this read
			auto read_it = r.refCoordToAP.find( pos );
			if ( read_it == r.refCoordToAP.end() ) continue;

			double obs = read_it->second->getSignalMean();
			if ( obs <= 0.0 or obs >= 250.0 ) continue;

			// Reference Gaussian at this position
			auto ref_it = contigRef.find( pos );
			if ( ref_it == contigRef.end() ) continue;

			const RefSigEntry &re = ref_it->second;
			if ( re.std <= 0.0 ) continue;    // single-read positions have std = 0

			double z = ( obs - re.mean ) / re.std;
			window_ll += -0.5 * ( LOG2PI + 2.0 * log( re.std ) + z * z );
			valid++;
		}

		if ( valid >= minWindowCov ){
			out += std::to_string( center ) + "\t"
			     + std::to_string( window_ll / static_cast<double>(valid) ) + "\t"
			     + std::to_string( valid ) + "\n";
		}
	}
}


int refScan_main( int argc, char** argv ){

	RefScanArguments args = parseRefScanArguments( argc, argv );

	//load DNAscent index
	std::map< std::string, IndexEntry > readID2path;
	parseIndex( args.indexFilename, readID2path );

	//load the reference signal profile
	std::cout << "Loading reference signal... ";
	std::map< std::string, std::map< int, RefSigEntry > > refSig = loadRefSignal( args.refSignalFilename );
	std::cout << "ok." << std::endl;

	//import fasta reference
	std::map< std::string, std::string > reference = import_reference_pfasta( args.referenceFilename );

	//open the bam
	std::cout << "Opening bam file... ";
	htsFile   *bam_fh_cr  = sam_open( args.bamFilename.c_str(), "r" );
	if (bam_fh_cr == NULL) throw IOerror( args.bamFilename );
	bam_hdr_t *bam_hdr_cr = sam_hdr_read( bam_fh_cr );
	std::cout << "ok." << std::endl;

	//open the output file
	std::ofstream outFile( args.outputFilename );
	if ( not outFile.is_open() ) throw IOerror( args.outputFilename );

	//write output header
	outFile << "#Alignment "      << args.bamFilename       << "\n";
	outFile << "#Genome "         << args.referenceFilename << "\n";
	outFile << "#Index "          << args.indexFilename     << "\n";
	outFile << "#RefSignal "      << args.refSignalFilename << "\n";
	outFile << "#Threads "        << args.threads           << "\n";
	outFile << "#MappingQuality " << args.minQ              << "\n";
	outFile << "#MappingLength "  << args.minL              << "\n";
	outFile << "#MinWindowCov "   << args.minWindowCov      << "\n";
	outFile << "#Version "        << VERSION                << "\n";
	outFile << "#Commit "         << getGitCommit()         << "\n";

	//open a log file
	std::cout << "Opening log file... ";
	std::string logFilename = strip_extension( args.outputFilename ) + ".refScan.log";
	std::ofstream logfile( logFilename );
	if (logfile.is_open()) std::cout << "ok." << std::endl;
	else throw IOerror( logFilename );

	//initialise progress
	int numOfRecords = 0, prog = 0, failed = 0;
	countRecords( bam_fh_cr, bam_hdr_cr, numOfRecords, args.minQ, args.minL );
	progressBar pb( numOfRecords, true );
	bam_hdr_destroy( bam_hdr_cr );
	hts_close( bam_fh_cr );

	std::mutex out_mtx;
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

			#pragma omp parallel for schedule(dynamic) shared(buffer,Pore_Substrate_Config,args,prog,failed,refSig,out_mtx,log_mtx) num_threads(args.threads)
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

				// Score the rolling k-mer windows for this read
				std::string readOut;
				scoreWindowsOnRead( r, refSig, args.minWindowCov, readOut );

				prog++;

				#pragma omp critical
				{
					outFile << readOut;
					pb.displayProgress( prog, failed, failedEvents );
				}
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
	outFile.close();
	logfile.close();
	return 0;
}
