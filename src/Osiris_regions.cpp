//----------------------------------------------------------
// Copyright 2017 University of Oxford
// Written by Michael A. Boemo (michael.boemo@path.ox.ac.uk)
// This software is licensed under GPL-3.0.  You should have
// received a copy of the license with this software.  If
// not, please Email the author.
//----------------------------------------------------------

#include <fstream>
#include "Osiris_regions.h"
#include "../Penthus/src/probability.h"
#include "../Penthus/src/hmm.h"
#include "common.h"
#include "data_IO.h"
#include "error_handling.h"
#include <cmath>

#define _USE_MATH_DEFINES

static const char *help=
"build: Osiris executable that builds detection data for later processing by Osiris detect.\n"
"To run Osiris regions, do:\n"
"  ./Osiris regions [arguments]\n"
"Example:\n"
"  ./Osiris regions -d /path/to/osiris_detect_output.out -r path/to/reference.fasta -o /path/to/output_prefix -t 20\n"
"Required arguments are:\n"
"  -d,--detect               path to output file from Osiris detect,\n"
"  -r,--reference            path to genome reference in fasta format,\n"
"  -o,--output               path to output bed prefix.\n"
"Optional arguments are:\n"
"  -t,--threads              number of threads (default is 1 thread).\n";

struct Arguments {
	std::string detectFilename;
	std::string outputFilename;
	std::string referenceFilename;
	int threads;
};

Arguments parseRegionsArguments( int argc, char** argv ){

	if( argc < 2 ){

		std::cout << "Exiting with error.  Insufficient arguments passed to Osiris build." << std::endl << help << std::endl;
		exit(EXIT_FAILURE);
	}

	if ( std::string( argv[ 1 ] ) == "-h" or std::string( argv[ 1 ] ) == "--help" ){

		std::cout << help << std::endl;
		exit(EXIT_SUCCESS);
	}
	else if( argc < 4 ){

		std::cout << "Exiting with error.  Insufficient arguments passed to Osiris regions." << std::endl;
		exit(EXIT_FAILURE);
	}

	Arguments args;

	/*defaults - we'll override these if the option was specified by the user */
	args.threads = 1;

	/*parse the command line arguments */
	for ( int i = 1; i < argc; ){

		std::string flag( argv[ i ] );

		if ( flag == "-d" or flag == "--detect" ){

			std::string strArg( argv[ i + 1 ] );
			args.detectFilename = strArg;
			i+=2;
		}
		else if ( flag == "-o" or flag == "--output" ){

			std::string strArg( argv[ i + 1 ] );
			args.outputFilename = strArg + ".bed";
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
		else throw InvalidOption( flag );
	}
	return args;
}


void callRegions( readDetection &rd ){

	HiddenMarkovModel hmm = HiddenMarkovModel( 2, 3 );

	SilentDistribution sd( 0.0, 0.0 );
	BernoulliDistribution distForB( 0.25, 0 );
	BernoulliDistribution distForT( 0.01, 0 );

	State BrdUregion( &distForB, "BrdU", "", "", 1.0);
	State Thymregion( &distForT, "Thym", "", "", 1.0);

	hmm.add_state(BrdUregion);
	hmm.add_state(Thymregion);

	/*from start */
	hmm.add_transition( hmm.start, BrdUregion, 0.5 );
	hmm.add_transition( hmm.start, Thymregion, 0.5 );

	/*from BrdU */
	hmm.add_transition( BrdUregion, BrdUregion, 1.0 - 2.0/100.0 );
	hmm.add_transition( BrdUregion, Thymregion, 1.0/100.0 );
	hmm.add_transition( BrdUregion, hmm.end, 1.0/100.0 );

	/*from Thym */
	hmm.add_transition( Thymregion, Thymregion, 1.0 - 2.0/100.0 );
	hmm.add_transition( Thymregion, BrdUregion, 1.0/100.0 );
	hmm.add_transition( Thymregion, hmm.end, 1.0/100.0 );

	hmm.finalise();

	std::pair< double, std::vector< std::string > > path = hmm.viterbi( rd.BrdUProb );

	//print out path for debugging
	//for ( int i = 0; i < path.second.size(); i++ ){

	//	std::cout << path.second[i] << " ";
	//}
	//std::cout << std::endl;

	int start, end;
	for ( int i = 1; i < path.second.size() - 1; i++ ){

		if ( path.second[i] == "BrdU" ){

			if (path.second[i-1] == "START" or path.second[i-1] == "Thym"){

				start = rd.positions[i];
			}
		}
		else{

			if (path.second[i-1] == "BrdU"){

				end = rd.positions[i-1];
				Track tr;
				//if ( end - start > 1000){
				tr.lowerBound = start;
				tr.upperBound = end;
				rd.tracks.push_back(tr);
				//}
			}
		}
	}
}


void writePSL( readDetection &rd, std::map< std::string, std::string > &reference){

	if (rd.tracks.size() == 0) return;

	std::cout << 0 << " "; //matches
	std::cout << 0 << " "; //mismatches
	std::cout << 0 << " "; //repMatches
	std::cout << 0 << " "; //nCount
	std::cout << 0 << " "; //qNumInsert 
	std::cout << 0 << " "; //qBaseInsert 
	std::cout << 0 << " "; //tNumInsert 
	std::cout << 0 << " "; //tBaseInsert 
	std::cout << "+" << " "; //strand
	std::cout << rd.readID << " "; //queryName
	std::cout << rd.mappingUpper - rd.mappingLower << " "; //qSize
	std::cout << 0 << " "; //qStart
	std::cout << rd.mappingUpper - rd.mappingLower << " "; //qEnd
	std::cout << rd.chromosome << " "; //tName
	std::cout << reference[rd.chromosome].size() << " "; //tSize
	std::cout << rd.mappingLower << " "; //tStart
	std::cout << rd.mappingUpper << " "; //tEnd
	std::cout << rd.tracks.size() << " "; //blockCount

	//blockSizes
	for ( int i = 0; i < rd.tracks.size(); i++ ){

		std::cout << rd.tracks[i].upperBound - rd.tracks[i].lowerBound  << ",";
	}
	std::cout << " ";

	//qStarts
	for ( int i = 0; i < rd.tracks.size(); i++ ){

		std::cout << rd.tracks[i].lowerBound - rd.mappingLower << ",";
	}
	std::cout << " ";

	//tStarts
	for ( int i = 0; i < rd.tracks.size(); i++ ){

		std::cout << rd.tracks[i].lowerBound << ",";
	}
	std::cout << std::endl;
}


int regions_main( int argc, char** argv ){

	Arguments args = parseRegionsArguments( argc, argv );

	std::map< std::string, std::string > reference = import_reference(args.referenceFilename);

	std::ifstream inFile( args.detectFilename );
	if ( not inFile.is_open() ) throw IOerror( args.detectFilename );

	std::ofstream outFile( args.outputFilename );
	if ( not outFile.is_open() ) throw IOerror( args.outputFilename );

	std::string line;
	std::vector< readDetection > buffer;
	int maxDistance = 500, starting = -1, gap = 0;
	while ( std::getline( inFile, line ) ){

		if ( line.substr(0,1) == ">" ){

			if ( buffer.size() > args.threads ){

				for ( int i = 0; i < buffer.size(); i++ ){

					if ( buffer[i].BrdUProb.size() > 0 ){

						callRegions( buffer[i] );
						writePSL( buffer[i], reference );
					}
				}
				buffer.clear();
			}

			readDetection rd;
			rd.readID = line.substr(1, line.find(' ') - 1);
			rd.chromosome = line.substr(line.find(' ') + 1, line.find(':') - line.find(' ') - 1);
			rd.mappingLower = std::stoi(line.substr(line.find(':') + 1, line.rfind('-') - line.find(':') - 1));
			rd.mappingUpper = std::stoi(line.substr(line.rfind('-')+1));
			buffer.push_back(rd);
		}
		else if ( distance > maxDistance ){




		}
		else{

			std::string column;
			std::stringstream ssLine(line);
			int position, cIndex = 0;
			while( std::getline( ssLine, column, '\t' ) ){

				if (cIndex == 0){

					if ( starting == -1 ){

						starting = std::stoi(column);
					}
					else{

						gap += std::stoi(column) - starting;
					}
				}
				else if (cIndex == 1){

					if ( std::stof(column) > 2.5 ) calls += 1;
					attempts += 1;



					buffer.back().BrdUProb.push_back( std::stoi( column ) );
				}
				cIndex++;
			}
		}
	}

	return 0;
}
