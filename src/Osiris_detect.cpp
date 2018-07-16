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
#include "../htslib/htslib/hts.h"
#include "../htslib/htslib/sam.h"
#include "../fast5/include/fast5.hpp"


static const char *help=
"build: Osiris executable that builds detection data for later processing by Osiris detect.\n"
"To run Osiris build, do:\n"
"  ./Osiris build [arguments]\n"
"Example:\n"
"  ./Osiris build -b /path/to/alignment.bam -r /path/to/reference.fasta -i /path/to/index.index -o /path/to/output.out -t 20\n"
"Required arguments are:\n"
"  -b,--bam                  path to alignment BAM file,\n"
"  -r,--reference            path to genome reference in fasta format,\n"
"  -m,--analogue-model       path to analogue model file,\n"
"  -i,--index                path to Osiris index,\n"
"  -o,--output               path to output file that will be generated.\n"
"Optional arguments are:\n"
"  -t,--threads              number of threads (default is 1 thread).\n";

struct Arguments {
	std::string bamFilename;
	std::string referenceFilename;
	std::string outputFilename;
	std::string indexFilename;
	std::string analogueModelFilename;
	int threads;
};

Arguments parseBuildArguments( int argc, char** argv ){

	if( argc < 2 ){

		std::cout << "Exiting with error.  Insufficient arguments passed to Osiris build." << std::endl << help << std::endl;
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

	Arguments args;

	/*defaults - we'll override these if the option was specified by the user */
	args.threads = 1;

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
		else if ( flag == "-m" or flag == "--analogue-model" ){

			std::string strArg( argv[ i + 1 ] );
			args.analogueModelFilename = strArg;
			i+=2;
		}
		else throw InvalidOption( flag );
	}
	return args;
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

			nd.push_back( NormalDistribution( analogueModel.at(sixMer).first, analogueModel.at(sixMer).second ) );
		}
		else {

			nd.push_back( NormalDistribution( SixMer_model.at(sixMer).first, SixMer_model.at(sixMer).second ) );
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



// from scrappie
float fast5_read_float_attribute(hid_t group, const char *attribute) {
    float val = NAN;
    if (group < 0) {
#ifdef DEBUG_FAST5_IO
        fprintf(stderr, "Invalid group passed to %s:%d.", __FILE__, __LINE__);
#endif
        return val;
    }

    hid_t attr = H5Aopen(group, attribute, H5P_DEFAULT);
    if (attr < 0) {
#ifdef DEBUG_FAST5_IO
        fprintf(stderr, "Failed to open attribute '%s' for reading.", attribute);
#endif
        return val;
    }

    H5Aread(attr, H5T_NATIVE_FLOAT, &val);
    H5Aclose(attr);

    return val;
}
//end scrappie


void parseCigar(bam1_t *record, std::map< int, int > &ref2query, int &refStart, int &refEnd ){

	//initialise reference and query coordinates for the first match
	refStart = record -> core.pos;
	int queryPosition = 0;
	int refPosition = 0;

	const uint32_t *cigar = bam_get_cigar(record);
	for ( int i = 0; i < record -> core.n_cigar; ++i){

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
	refEnd = refStart + refPosition;
}


void parseIndex( std::string indexFilename, std::map< std::string, std::string > &readID2path ){

	std::ifstream indexFile( indexFilename );
	if ( not indexFile.is_open() ) throw IOerror( indexFilename );

	std::string line;
	while ( std::getline( indexFile, line) ){

		std::string readID = line.substr(0, line.find('\t'));
		std::string path = line.substr(line.find('\t')+1);
		readID2path[readID] = path;
	}
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

void getEvents( std::string fast5Filename, std::vector<double> &raw ){

	//open the file
	hid_t hdf5_file = H5Fopen(fast5Filename.c_str(), H5F_ACC_RDONLY, H5P_DEFAULT);

	//get the channel parameters
	const char *scaling_path = "/UniqueGlobalKey/channel_id";
	hid_t scaling_group = H5Gopen(hdf5_file, scaling_path, H5P_DEFAULT);
	float digitisation = fast5_read_float_attribute(scaling_group, "digitisation");
	float offset = fast5_read_float_attribute(scaling_group, "offset");
	float range = fast5_read_float_attribute(scaling_group, "range");
	float sample_rate = fast5_read_float_attribute(scaling_group, "sampling_rate");

	//get the raw signal
	hid_t space;
	hsize_t nsample;
	herr_t status;
	float raw_unit;
	float *rawptr = NULL;

	ssize_t size = H5Lget_name_by_idx(hdf5_file, "/Raw/Reads/", H5_INDEX_NAME, H5_ITER_INC, 0, NULL, 0, H5P_DEFAULT);
	char* name = (char*)calloc(1 + size, sizeof(char));
	H5Lget_name_by_idx(hdf5_file, "/Raw/Reads/", H5_INDEX_NAME, H5_ITER_INC, 0, name, 1 + size, H5P_DEFAULT);
	std::string readName(name);
	free(name);
	std::string signal_path = "/Raw/Reads/" + readName + "/Signal";

	hid_t dset = H5Dopen(hdf5_file, signal_path.c_str(), H5P_DEFAULT);
	if (dset < 0 ) throw BadFast5Field(); 
	space = H5Dget_space(dset);
	if (space < 0 ) throw BadFast5Field(); 
	H5Sget_simple_extent_dims(space, &nsample, NULL);
   	rawptr = (float*)calloc(nsample, sizeof(float));
    	status = H5Dread(dset, H5T_NATIVE_FLOAT, H5S_ALL, H5S_ALL, H5P_DEFAULT, rawptr);
	
	raw_unit = range / digitisation;
	for ( size_t i = 0; i < nsample; i++ ){

		raw.push_back( (rawptr[i] + offset) * raw_unit );
	}

	H5Fclose(hdf5_file);
}


int countRecords( htsFile *bam_fh, hts_idx_t *bam_idx, bam_hdr_t *bam_hdr, int &numOfRecords, double &avgRecordLength ){

	hts_itr_t* itr = sam_itr_querys(bam_idx,bam_hdr,".");
	int result;
	int recordLengths = 0;
	do {

		bam1_t *record = bam_init1();
		result = sam_itr_next(bam_fh, itr, record);
		numOfRecords++;
		recordLengths += record -> core.l_qseq;
	} while (result > 0);
	avgRecordLength = (double)recordLengths / (double)numOfRecords;
}


std::vector< unsigned int > getPOIs( std::string &refSeq, std::map< std::string, std::pair< double, double > > &analogueModel, int windowLength ){

	std::vector< unsigned int > POIs;

	for ( unsigned int i = windowLength; i < refSeq.length() - 2*windowLength; i++ ){

		if ( analogueModel.count( refSeq.substr(i, 6) ) > 0 ) POIs.push_back(i);
	}
	return POIs;
}


int detect_main( int argc, char** argv ){

	Arguments args = parseBuildArguments( argc, argv );

	//load osiris index
	std::map< std::string, std::string > readID2path;
	parseIndex( args.indexFilename, readID2path );

	//import fasta reference
	std::map< std::string, std::string > reference = import_reference( args.referenceFilename );

	/*import the analogue pore model that we specified on the command line */
	std::map< std::string, std::pair< double, double > > analogueModel =  import_poreModel( args.analogueModelFilename );

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
	double avgRecordLength;
	countRecords( bam_fh, bam_idx, bam_hdr, numOfRecords, avgRecordLength );
	progressBar pb(numOfRecords);

	//build an iterator for all reads in the bam file
	const char *allReads = ".";
	itr = sam_itr_querys(bam_idx,bam_hdr,allReads);

	int windowLength = 20;
	int result;
	do {

		read r;

		//initialise the record and get the record from the file iterator
		bam1_t *record = bam_init1();
		result = sam_itr_next(bam_fh, itr, record);

		//skip reverse complements
		if ( bam_is_rev(record) ) continue;
		
		//iterate on the cigar string to fill up the reference-to-query coordinate map
		int refStart, refEnd;
		parseCigar(record, r.refToQuery, refStart, refEnd);

		//fetch the basecall from the bam file
		r.basecall = getQuerySequence(record);

		//get the read name (which will be the ONT readID from Albacore basecall)
		const char *queryName = bam_get_qname(record);
		std::string s_queryName(queryName);
		r.readID = s_queryName;

		//get the name of the reference mapped to
		std::string mappedTo(bam_hdr -> target_name[record -> core.tid]);
		r.referenceMappedTo = mappedTo;

		//open fast5 and normalise events to pA
		r.filename = readID2path[s_queryName];
		try{
			
			getEvents( r.filename, r.raw );
		}
		catch ( BadFast5Field &bf5 ){
			
			failed++;
			prog++;
			continue;
		}
		normaliseEvents(r);

		/*disregard this event if the quality score is too low */
		if ( fabs(r.qualityScore) > 2 ){

			failed++;
			prog++;
			continue;
		}

		/*get the subsequence of the reference this read mapped to and build an HMM from it */
		std::string refSeqMapped = reference[r.referenceMappedTo].substr(refStart, refEnd-refStart);
			
		/*push the filename for this read to the output */
		std::stringstream ss;
		ss << ">" << r.readID << " " << r.referenceMappedTo << ":" << refStart << "-" << refEnd << std::endl;

		std::vector< unsigned int > POIs = getPOIs( refSeqMapped, analogueModel, windowLength );
		#pragma omp parallel for ordered schedule(static) shared(ss,POIs,failed,pb,prog,SixMer_model,analogueModel) num_threads(args.threads)
		for ( unsigned int i = 0; i < POIs.size(); i++ ){

			int BrdULoc = POIs[i];
			int posOnQuery = (r.refToQuery).at(BrdULoc);

			std::string readSnippet = refSeqMapped.substr(BrdULoc - windowLength, 2*windowLength);
			std::vector< double > eventSnippet;

			/*get the events that correspond to the read snippet */
			for ( unsigned int j = 0; j < (r.eventAlignment).size(); j++ ){

				/*if an event has been aligned to a position in the window, add it */
				if ( (r.eventAlignment)[j].second >= (r.refToQuery)[BrdULoc - windowLength] and (r.eventAlignment)[j].second <= (r.refToQuery)[BrdULoc + windowLength] ){

					eventSnippet.push_back( (r.normalisedEvents)[j] );
				}

				/*stop once we get to the end of the window */
				if ( (r.eventAlignment)[j].second > (r.refToQuery)[BrdULoc + windowLength] ) break;
			}
			double logProbThymidine = seqProbability(readSnippet, eventSnippet, analogueModel, windowLength, false);
			double logProbAnalogue = seqProbability(readSnippet, eventSnippet, analogueModel, windowLength, true);
			double logLikelihoodRatio = logProbAnalogue - logProbThymidine;

			#pragma omp ordered
			ss << BrdULoc + refStart << "\t" << logLikelihoodRatio << "\t" << refSeqMapped.substr(BrdULoc, 6) << "\t" << (r.basecall).substr(posOnQuery, 6) << std::endl;
		}
		outFile << ss.rdbuf();
		prog++;
		pb.displayProgress( prog, failed );
	
	} while (result > 0);

	return 0;
}
