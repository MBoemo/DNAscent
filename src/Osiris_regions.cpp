//----------------------------------------------------------
// Copyright 2017 University of Oxford
// Written by Michael A. Boemo (michael.boemo@path.ox.ac.uk)
// This software is licensed under GPL-3.0.  You should have
// received a copy of the license with this software.  If
// not, please Email the author.
//----------------------------------------------------------
 #include <fstream>
#include "Osiris_regions.h"
#include "data_IO.h"
#include "error_handling.h"
#include <cmath>
#include <math.h>
 #define _USE_MATH_DEFINES

 static const char *help=
"build: Osiris executable that builds a PSL file from the output of Osiris detect.\n"
"To run Osiris regions, do:\n"
"  ./Osiris regions [arguments]\n"
"Example:\n"
"  ./Osiris regions -d /path/to/osiris_detect_output.out -r path/to/reference.fasta -o /path/to/output_prefix -t 20\n"
"Required arguments are:\n"
"  -d,--detect               path to output file from Osiris detect,\n"
"  -p,--probability          probability that a thymidine 6mer contains a BrdU,\n"
"  -o,--output               path to output directory for bedgraph files.\n"
"Optional arguments are:\n"
"  -r,--resolution           minimum length of regions (default is 2kb).\n"
"  -z,--zScore               zScore threshold for BrdU call (default is -2).\n";

 struct Arguments {

	std::string detectFilename;
	double probability, threshold;
	unsigned int resolution;
	std::string outputFilename;
};

 Arguments parseRegionsArguments( int argc, char** argv ){

 	if( argc < 2 ){
 		std::cout << "Exiting with error.  Insufficient arguments passed to Osiris regions." << std::endl << help << std::endl;
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
	args.resolution = 2000;
	args.threshold = -2;

 	/*parse the command line arguments */
	for ( int i = 1; i < argc; ){

 		std::string flag( argv[ i ] );

 		if ( flag == "-d" or flag == "--detect" ){
 			std::string strArg( argv[ i + 1 ] );
			args.detectFilename = strArg;
			i+=2;
		}
		else if ( flag == "-p" or flag == "--probability" ){
 			std::string strArg( argv[ i + 1 ] );
			args.probability = std::stof( strArg.c_str() );
			i+=2;
		}
		else if ( flag == "-o" or flag == "--output" ){
 			std::string strArg( argv[ i + 1 ] );
			args.outputFilename = strArg;
			i+=2;
		}
		else if ( flag == "-z" or flag == "--zScore" ){
 			std::string strArg( argv[ i + 1 ] );
			args.threshold = std::stof(strArg.c_str());
			i+=2;
		}
		else if ( flag == "-r" or flag == "--resolution" ){
 			std::string strArg( argv[ i + 1 ] );
			args.resolution = std::stoi( strArg.c_str() );
			i+=2;
		}
		else throw InvalidOption( flag );
	}
	return args;
}


struct region{

	std::string call;
	int start, end;
	double score;
};


int regions_main( int argc, char** argv ){

	Arguments args = parseRegionsArguments( argc, argv );

 	std::ifstream inFile( args.detectFilename );
	if ( not inFile.is_open() ) throw IOerror( args.detectFilename );

 	std::ofstream outFile( args.outputFilename );
	if ( not outFile.is_open() ) throw IOerror( args.outputFilename );

 	std::string line, header;
	std::vector< region > buffer;
	int calls = 0, attempts = 0, gap = 0, startingPos = -1;
	while( std::getline( inFile, line ) ){

		if ( line.substr(0,1) == ">" ){

			if ( buffer.size() > 0 ){

				outFile << header << std::endl;
				
				for ( auto r = buffer.begin(); r < buffer.end(); r++ ){

					outFile << r -> start << "\t" << r -> end << "\t" << r -> score << "\t" << r -> call <<  std::endl;
				}
			}
			header = line;
			buffer.clear();
			calls = 0, attempts = 0, gap = 0, startingPos = -1;
			
		}
		else{

			std::string column;
			std::stringstream ssLine(line);
			int position, cIndex = 0;
			double B;
			while ( std::getline( ssLine, column, '\t' ) ){

				if ( cIndex == 0 ){

					position = std::stoi(column);
				}
				else if ( cIndex == 1 ){

					B = std::stof(column);
					if ( B > 2.5 ){

						calls++;
					}
				}
				cIndex++;
			}

			if ( startingPos == -1 ) startingPos = position;
			gap = position - startingPos;
			attempts++;

			if ( gap > args.resolution ){

				region r;

				r.score = (calls - attempts * args.probability) / sqrt( attempts * args.probability * ( 1 - args.probability) );

				if ( r.score > args.threshold ) r.call = "BrdU";
				else r.call = "Thym";

				r.start = startingPos;
				r.end = position;
				
				buffer.push_back(r);
				calls = 0, attempts = 0, gap = 0, startingPos = -1;
			}
		}
	}
	return 0;
}
