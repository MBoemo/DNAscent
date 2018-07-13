//----------------------------------------------------------
// Copyright 2017 University of Oxford
// Written by Michael A. Boemo (michael.boemo@path.ox.ac.uk)
// This software is licensed under GPL-3.0.  You should have
// received a copy of the license with this software.  If
// not, please Email the author.
//----------------------------------------------------------

#include <fstream>
#include "Osiris_build.h"
#include "common.h"
#include "build_model.h"
#include "data_IO.h"
#include "error_handling.h"
#include "event_handling.h"
#include "poreModels.h"
#include "poreSpecificParameters.h"
#include "../htslib-1.8/htslib/hts.h"
#include "../htslib-1.8/htslib/sam.h"
#include "../fast5/include/fast5.hpp"


static const char *help=
"build: Osiris executable that builds detection data for later processing by Osiris detect.\n"
"To run Osiris build, do:\n"
"  ./Osiris build [arguments]\n"
"Example:\n"
"  ./Osiris build -b /path/to/alignment.bam -r /path/to/reference.fasta -i /path/to/index.index -o outputPrefix -t 20\n"
"Required arguments are:\n"
"  -b,--bam                  path to alignment BAM file,\n"
"  -r,--reference            path to genome reference in fasta format,\n"
"  -i,--index                path to Osiris index.\n"
"Optional arguments are:\n"
"  -t,--threads              number of threads (default is 1 thread).\n";

struct Arguments {
	std::string bamFilename;
	std::string referenceFilename;
	std::string outputFilename;
	std::string indexFilename;
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
		else throw InvalidOption( flag );
	}
	return args;
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


void parseCigar(bam1_t *record, std::map< int, int > &ref2query ){

	//initialise reference and query coordinates for the first match
	uint32_t refPosition = record -> core.pos;
	int queryPosition = 0;

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
		//for a deletion, advance only the reference position (N.B. we're not worried about insertions)
		else if (op == BAM_CDEL or op == BAM_CREF_SKIP){

			for ( int j = refPosition; j < refPosition + ol; j++ ){

				ref2query[j] = queryPosition;
			}
			refPosition += ol;
		}
	}
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

	std::cout << digitisation << std::endl;

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
	//if dset < 0 catch 
	space = H5Dget_space(dset);
	//if space < 0 catch
	H5Sget_simple_extent_dims(space, &nsample, NULL);
   	rawptr = (float*)calloc(nsample, sizeof(float));
    	status = H5Dread(dset, H5T_NATIVE_FLOAT, H5S_ALL, H5S_ALL, H5P_DEFAULT, rawptr);
	
	raw_unit = range / digitisation;
	for ( size_t i = 0; i < nsample; i++ ){

		raw.push_back( (rawptr[i] + offset) * raw_unit );
	}

	H5Fclose(hdf5_file);
}


int build_main( int argc, char** argv ){

	Arguments args = parseBuildArguments( argc, argv );

	std::map< std::string, std::string > readID2path;
	parseIndex( args.indexFilename, readID2path );

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

	//build an iterator for all reads in the bam file
	const char *allReads = ".";
	itr = sam_itr_querys(bam_idx,bam_hdr,allReads);

	int result;
	do {

		read r;

		//initialise the record and get the record from the file iterator
		bam1_t *record = bam_init1();
		result = sam_itr_next(bam_fh, itr, record);
		
		//iterate on the cigar string to fill up the reference-to-query coordinate map
		parseCigar(record, r.refToQuery);

		//fetch the basecall from the bam file
		r.basecall = getQuerySequence(record);

		//get the read name (which will be the ONT readID from Albacore basecall)
		const char *queryName = bam_get_qname(record);
		std::string s_queryName(queryName);
		r.readID = s_queryName;

		//open fast5 and normalise events to pA
		float *rawSignal = NULL;
		getEvents( readID2path[s_queryName], r.raw );
		normaliseEvents(r);

		//run HMM detection on this read
		
	} while (result > 0);

	return 0;
}
