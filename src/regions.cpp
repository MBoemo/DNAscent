//----------------------------------------------------------
// Copyright 2017 University of Oxford
// Written by Michael A. Boemo (michael.boemo@path.ox.ac.uk)
// This software is licensed under GPL-2.0.  You should have
// received a copy of the license with this software.  If
// not, please Email the author.
//----------------------------------------------------------

 #include <fstream>
#include "regions.h"
#include "data_IO.h"
#include "error_handling.h"
#include <cmath>
#include <math.h>
 #define _USE_MATH_DEFINES

 static const char *help=
"build: DNAscent executable that builds a PSL file from the output of DNAscent detect.\n"
"To run DNAscent regions, do:\n"
"  ./DNAscent regions [arguments]\n"
"Example:\n"
"  ./DNAscent regions -d /path/to/regions_output.out -p 0.2 -o /path/to/output_prefix\n"
"Required arguments are:\n"
"  -d,--detect               path to output file from DNAscent detect,\n"
"  -p,--probability          probability that a thymidine 6mer contains a BrdU,\n"
"  -o,--output               path to output directory for bedgraph files.\n"
"Optional arguments are:\n"
"     --replication          detect fork direction and call origin firing (default: off),\n"
"  -r,--resolution           minimum length of regions (default is 2kb).\n"
"  -z,--zScore               zScore threshold for BrdU call (default is -2).\n";

 struct Arguments {

	std::string detectFilename;
	double probability, threshold;
	unsigned int resolution;
	std::string outputFilename;
	bool callReplication;
};

Arguments parseRegionsArguments( int argc, char** argv ){

 	if( argc < 2 ){
 		std::cout << "Exiting with error.  Insufficient arguments passed to DNAscent regions." << std::endl << help << std::endl;
		exit(EXIT_FAILURE);
	}
 	if ( std::string( argv[ 1 ] ) == "-h" or std::string( argv[ 1 ] ) == "--help" ){
 		std::cout << help << std::endl;
		exit(EXIT_SUCCESS);
	}
	else if( argc < 4 ){
 		std::cout << "Exiting with error.  Insufficient arguments passed to DNAscent regions." << std::endl;
		exit(EXIT_FAILURE);
	}

 	Arguments args;

 	/*defaults - we'll override these if the option was specified by the user */
	args.resolution = 2000;
	args.threshold = -2;
	args.callReplication = false;

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
		else if ( flag == "--replication" ){

			args.callReplication = true;
			i+=1;
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
	std::string forkDir;
};


void callOrigins( std::vector< region > &regions, double threshold, std::string &header, std::ofstream &repFile ){

	//10kb moving average filter
	std::vector< double > out;
	for ( int i = 2; i < regions.size() - 2; i++ ){

		out.push_back( ( regions[i-2].score + regions[i-1].score + regions[i].score + regions[i+1].score + regions[i+2].score ) / 5.0 );
	}
	for ( int i = 0; i < out.size(); i++ ){

		regions[i+2].score = out[i];
	}

	//secondary 6kb moving average filter
	out.clear();
	for ( int i = 1; i < regions.size() - 1; i++ ){

		out.push_back( ( regions[i-1].score + regions[i].score + regions[i+1].score ) / 3.0 );
	}
	for ( int i = 0; i < out.size(); i++ ){

		regions[i+1].score = out[i];
		
		//re-assess the call with the smoothed scores
		if ( regions[i+1].score >= 0 ) regions[i+1].call = "BrdU";
		else regions[i+1].call = "Thym";
	}
	out.clear();

	std::vector< double > derivatives;
	
	//get the central derivative for all the middle regions
	for ( int i = 1; i < regions.size() - 1; i++ ){

		double former = regions[i-1].score;
		double next = regions[i+1].score;

		derivatives.push_back( ( next - former ) / 2.0 );
	}

	//error correct for derivatives that are sandwiched between two of the opposite sign
	for ( int i = 2; i < derivatives.size() - 2; i++ ){

		if ( derivatives[i-1] < 0 and derivatives[i+1] < 0 and derivatives[i] > 0 ) derivatives[i] = ( derivatives[i-1] + derivatives[i+1] ) / 2.0;
		if ( derivatives[i-1] > 0 and derivatives[i+1] > 0 and derivatives[i] < 0 ) derivatives[i] = ( derivatives[i-1] + derivatives[i+1] ) / 2.0;
	}

	//smooth derivative calls
	for ( int i = 1; i < derivatives.size() - 1; i++ ){

		out.push_back( ( derivatives[i-1] + derivatives[i] + derivatives[i+1] ) / 3.0 );
	}
	for ( int i = 0; i < out.size(); i++ ){

		derivatives[i+1] = out[i];
	}
	
	//error correct on BrdU calls for the middle
	for ( int i = 1; i < regions.size() - 1; i++ ){

		if ( regions[i-1].call == "BrdU" and regions[i+1].call == "BrdU" and regions[i].call == "Thym" ) regions[i].call = "BrdU";
		if ( regions[i-1].call == "Thym" and regions[i+1].call == "Thym" and regions[i].call == "BrdU" ) regions[i].call = "Thym";

	}

	//error correct on BrdU calls at the ends
	if ( regions.size() > 2 ){
	
		if ( regions[0].call == "BrdU" and regions[1].call == "Thym" and regions[2].call == "Thym" ) regions[0].call = "Thym";
		if ( regions[0].call == "Thym" and regions[1].call == "BrdU" and regions[2].call == "BrdU" ) regions[0].call = "BrdU";

		if ( regions[regions.size()-1].call == "BrdU" and regions[regions.size()-2].call == "Thym" and regions[regions.size()-3].call == "Thym" ) regions[0].call = "Thym";
		if ( regions[regions.size()-1].call == "Thym" and regions[regions.size()-2].call == "BrdU" and regions[regions.size()-3].call == "BrdU" ) regions[0].call = "BrdU";
	}

	//set the fork direction
	for ( int i = 1; i < regions.size() - 1; i++ ){

		if ( regions[i].call == "BrdU" and derivatives[i] < 0 ) regions[i].forkDir = "+";
		if ( regions[i].call == "BrdU" and derivatives[i] > 0 ) regions[i].forkDir = "-";
	}

	//call origin positions
	std::vector< int > oriPosForRead;
	for ( int i = 1; i < regions.size()-2; i++ ){

		if ( regions[i-1].forkDir == "-" and regions[i].forkDir == "-" and regions[i+1].forkDir == "+" and regions[i+2].forkDir == "+" ){

			if ( regions[i].score > regions[i+1].score ) oriPosForRead.push_back( (regions[i].start + regions[i].end) / 2 );
			else oriPosForRead.push_back( (regions[i+1].start + regions[i+1].end) / 2 );
		}
	}

	//write the origins to the output stream
	if ( oriPosForRead.size() > 0 ){

		repFile << header << std::endl;

		for ( unsigned int i = 0; i < oriPosForRead.size(); i++ ){		
		
			repFile << oriPosForRead[i] << std::endl;
		}
	}
}


int regions_main( int argc, char** argv ){

	Arguments args = parseRegionsArguments( argc, argv );

 	std::ifstream inFile( args.detectFilename );
	if ( not inFile.is_open() ) throw IOerror( args.detectFilename );

 	std::ofstream outFile( args.outputFilename );
	if ( not outFile.is_open() ) throw IOerror( args.outputFilename );

	std::ofstream repFile;
	if ( args.callReplication ){

		repFile.open("calledOrigins.dnascent");
	}


 	std::string line, header;
	std::vector< region > buffer;
	int calls = 0, attempts = 0, gap = 0, startingPos = -1;
	while( std::getline( inFile, line ) ){

		if ( line.substr(0,1) == ">" ){

			if ( buffer.size() > 5 ){

				outFile << header << std::endl;
			
				if (args.callReplication) callOrigins(buffer,args.threshold, header, repFile);
				
				for ( auto r = buffer.begin(); r < buffer.end(); r++ ){

					outFile << r -> start << "\t" << r -> end << "\t" << r -> score << "\t" << r -> call << "\t" << r -> forkDir <<  std::endl;
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

			if ( gap > args.resolution and attempts >= args.resolution / 30 ){

				region r;

				r.score = (calls - attempts * args.probability) / sqrt( attempts * args.probability * ( 1 - args.probability) );

				if ( r.score > args.threshold ) r.call = "BrdU";
				else r.call = "Thym";
				
				r.score += fabs(args.threshold);

				r.start = startingPos;
				r.end = position;
				
				buffer.push_back(r);
				calls = 0, attempts = 0, gap = 0, startingPos = -1;
			}
		}
	}

	//empty the buffer at the end
	if ( buffer.size() > 5 ){

		outFile << header << std::endl;
			
		if (args.callReplication) callOrigins(buffer,args.threshold, header, repFile);
				
		for ( auto r = buffer.begin(); r < buffer.end(); r++ ){

			outFile << r -> start << "\t" << r -> end << "\t" << r -> score << "\t" << r -> call << "\t" << r -> forkDir <<  std::endl;
		}
	}

	if ( repFile.is_open() ) repFile.close();
	inFile.close();
	outFile.close();

	return 0;
}