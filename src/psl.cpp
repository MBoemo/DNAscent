//----------------------------------------------------------
// Copyright 2017 University of Oxford
// Written by Michael A. Boemo (mb915@cam.ac.uk)
// This software is licensed under GPL-2.0.  You should have
// received a copy of the license with this software.  If
// not, please Email the author.
//----------------------------------------------------------

#include <fstream>
#include "psl.h"
#include "common.h"
#include "data_IO.h"
#include "error_handling.h"
#include <cmath>
#define _USE_MATH_DEFINES

 static const char *help=
"psl: DNAscent executable that builds a PSL file from the output of DNAscent detect.\n"
"To run DNAscent psl, do:\n"
"  ./DNAscent psl -d /path/to/DNAscentOutput.detect -r /path/to/reference.fasta -o /path/to/psl_prefix\n"
"Required arguments are:\n"
"  -d,--detect               path to output file from DNAscent detect,\n"
"  -r,--reference            path to genome reference in fasta format,\n"
"  -o,--output               path to output bed prefix.\n"
"Optional arguments are:\n"
"  -l,--likelihood           log-likelihood threshold for a positive analogue call (default: 1.25),\n"
"  -c,--cooldown             minimum gap between positive analogue calls (default: 4),\n"
"     --min                  minimum read length to compute (default is 1),\n"
"     --max                  maximum read length to compute (default is Inf).\n"
"Written by Michael Boemo, Department of Pathology, University of Cambridge.\n"
"Please submit bug reports to GitHub Issues (https://github.com/MBoemo/DNAscent/issues).";


 struct Arguments {
	std::string detectFilename;
	std::string outputFilename;
	std::string referenceFilename;
	bool cropToMin = false;
	unsigned int min = 0;
	bool cropToMax = false;
	unsigned int max = 0;
	double likelihood;
	int cooldown;


};
 Arguments parsePslArguments( int argc, char** argv ){
 	if( argc < 2 ){
 		std::cout << "Exiting with error.  Insufficient arguments passed to DNAscent psl." << std::endl << help << std::endl;
		exit(EXIT_FAILURE);
	}
 	if ( std::string( argv[ 1 ] ) == "-h" or std::string( argv[ 1 ] ) == "--help" ){
 		std::cout << help << std::endl;
		exit(EXIT_SUCCESS);
	}
	else if( argc < 4 ){
 		std::cout << "Exiting with error.  Insufficient arguments passed to DNAscent psl." << std::endl;
		exit(EXIT_FAILURE);
	}
 	Arguments args;
	args.likelihood = 1.25;
	args.cooldown = 4;

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
			args.outputFilename = strArg + ".psl";
			i+=2;
		}
		else if ( flag == "--min" ){
			args.cropToMin = true;
 			std::string strArg( argv[ i + 1 ] );
			args.min = std::stoi( strArg.c_str() );
			i+=2;
		}
		else if ( flag == "--max" ){
			args.cropToMax = true;
 			std::string strArg( argv[ i + 1 ] );
			args.max = std::stoi( strArg.c_str() );
			i+=2;
		}
		else if ( flag == "-r" or flag == "--reference" ){
 			std::string strArg( argv[ i + 1 ] );
			args.referenceFilename = strArg;
			i+=2;
		}
		else if ( flag == "-c" or flag == "--cooldown" ){
 			std::string strArg( argv[ i + 1 ] );
			args.cooldown = std::stoi( strArg.c_str() );
			i+=2;
		}
		else if ( flag == "-l" or flag == "--likelihood" ){
 			std::string strArg( argv[ i + 1 ] );
			args.likelihood = std::stof( strArg.c_str() );
			i+=2;
		}
		else throw InvalidOption( flag );
	}
	if (args.outputFilename == args.referenceFilename or args.outputFilename == args.detectFilename) throw OverwriteFailure();

	return args;
}

 void writePSL( readDetection &rd, std::map< std::string, std::string > &reference, std::ofstream &outFile ){
	if (rd.positions.size() == 0) return;
 	outFile << 0 << " "; //matches
	outFile << 0 << " "; //mismatches
	outFile << 0 << " "; //repMatches
	outFile << 0 << " "; //nCount
	outFile << 0 << " "; //qNumInsert 
	outFile << 0 << " "; //qBaseInsert 
	outFile << 0 << " "; //tNumInsert 
	outFile << 0 << " "; //tBaseInsert 
	outFile << rd.direction << " "; //strand
	outFile << rd.readID << " "; //queryName
	outFile << rd.mappingUpper - rd.mappingLower << " "; //qSize
	outFile << 0 << " "; //qStart
	outFile << rd.mappingUpper - rd.mappingLower << " "; //qEnd
	outFile << rd.chromosome << " "; //tName
	outFile << reference[rd.chromosome].size() << " "; //tSize
	outFile << rd.mappingLower << " "; //tStart
	outFile << rd.mappingUpper << " "; //tEnd
	outFile << rd.positions.size() + 2 << " "; //blockCount
 	//blockSizes
	outFile << 1 << ","; //extra for start
	for ( unsigned int i = 0; i < rd.positions.size(); i++ ){
 		outFile << 1  << ",";
	}
	outFile << 1 << ","; //extra for end
	outFile << " ";
 	//qStarts
	outFile << 0 << ","; //extra for start
	for ( unsigned int i = 0; i < rd.positions.size(); i++ ){
 		outFile << rd.positions[i] - rd.mappingLower << ",";
	}
	outFile << rd.mappingUpper - rd.mappingLower << ","; //extra for end
	outFile << " ";
 	//tStarts
	outFile << rd.mappingLower << ","; //extra for start
	for ( unsigned int i = 0; i < rd.positions.size(); i++ ){
 		outFile << rd.positions[i] << ",";
	}
	outFile << rd.mappingUpper << ","; //extra for end
	outFile << std::endl;
}
 int psl_main( int argc, char** argv ){

 	Arguments args = parsePslArguments( argc, argv );
 	std::map< std::string, std::string > reference = import_reference_pfasta(args.referenceFilename);

 	std::ifstream inFile( args.detectFilename );
	if ( not inFile.is_open() ) throw IOerror( args.detectFilename );
 	std::ofstream outFile( args.outputFilename );
	if ( not outFile.is_open() ) throw IOerror( args.outputFilename );

 	std::string line;
	std::vector< readDetection > buffer;
	bool recordRead = true;
	int callCooldown = 0;

	while ( std::getline( inFile, line ) ){

 		if ( line.substr(0,1) == ">" ){

 			if ( buffer.size() >= 10 ){

 				for ( unsigned int i = 0; i < buffer.size(); i++ ){

					writePSL( buffer[i], reference, outFile );
				}
				buffer.clear();
			}
 			readDetection rd;

			callCooldown = 0;

			std::stringstream ssLine(line);
			std::string column;
			int cIndex = 0;
			while ( std::getline( ssLine, column, ' ' ) ){

				if ( cIndex == 0 ) rd.readID = column;
				else if ( cIndex == 1 ) rd.chromosome = column;
				else if ( cIndex == 2 ) rd.mappingLower = std::stoi(column);
				else if ( cIndex == 3 ) rd.mappingUpper = std::stoi(column);
				else if ( cIndex == 4 ) rd.strand = column;
				else throw DetectParsing();
				cIndex++;
			}
			assert(rd.mappingUpper > rd.mappingLower);
			unsigned int readLength = rd.mappingUpper - rd.mappingLower;
			recordRead = true;
			if ( (args.cropToMin and readLength < args.min) or (args.cropToMax and readLength > args.max) ){

				recordRead = false;
				continue;
			}

			if (rd.strand == "fwd") rd.direction = "+";
			else if (rd.strand == "rev") rd.direction = "-";
			else throw BadStrandDirection();
			buffer.push_back(rd);
		}
		else if ( recordRead ){

			std::string column;
			std::stringstream ssLine(line);
			int position, cIndex = 0;
			double B;

			while( std::getline( ssLine, column, '\t' ) ){

				if (cIndex == 0){

					position = std::stoi(column);
				}
				else if (cIndex == 1){

					B = std::stof(column);

					if ( B > args.likelihood and position - callCooldown >= args.cooldown ){

						buffer.back().positions.push_back(position);
						callCooldown = position;
					}
					break;
				}
				cIndex++;
			}
		}
	}

 	for ( unsigned int i = 0; i < buffer.size(); i++ ){

		writePSL( buffer[i], reference, outFile );
	}
	buffer.clear();

 	return 0;
}
