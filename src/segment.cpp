//----------------------------------------------------------
// Copyright 2020 University of Cambridge
// Written by Michael A. Boemo (mb915@cam.ac.uk)
// This software is licensed under GPL-3.0.  You should have
// received a copy of the license with this software.  If
// not, please Email the author.
//----------------------------------------------------------

#include "../tensorflow/include/tensorflow/c/eager/c_api.h"
#include <fstream>
#include "data_IO.h"
#include "trainGMM.h"
#include "forkSense.h"
#include "segment.h"
#include "tensor.h"
#include <cmath>
#include <memory>
#include <math.h>
#include <algorithm>
#include <limits>
#include <stdlib.h>

static const char *help=
"segment: DNAscent executable that calls segments of high BrdU and EdU incorporation on each read.\n"
"DNAscent segment is designed to be a more general alternative to forkSense for users detecting\n"
"base analogues for some purpose other than DNA replication.\n"
"To run DNAscent segment, do:\n"
"   DNAscent segment -d /path/to/BrdU_EdU_calls.detect\n"
"Required arguments are:\n"
"  -d,--detect               path to output file from DNAscent detect,\n"
"Optional arguments are:\n"
"  -m,--minLen               length (in bp) of the shortest expected segment (default: 1000 bp),\n"
"  -t,--threads              number of threads (default: 1 thread).\n"
"Output is written to two bed files in the present working directory.\n"
"Written by Michael Boemo, Department of Pathology, University of Cambridge.\n"
"Please submit bug reports to GitHub Issues (https://github.com/MBoemo/DNAscent/issues).";


segmentArgs parseSegmentArguments( int argc, char** argv ){

 	if( argc < 2 ){
 		std::cout << "Exiting with error.  Insufficient arguments passed to DNAscent segment." << std::endl << help << std::endl;
		exit(EXIT_FAILURE);
	}
 	if ( std::string( argv[ 1 ] ) == "-h" or std::string( argv[ 1 ] ) == "--help" ){
 		std::cout << help << std::endl;
		exit(EXIT_SUCCESS);
	}
	else if( argc < 4 ){
 		std::cout << "Exiting with error.  Insufficient arguments passed to DNAscent segment." << std::endl;
		exit(EXIT_FAILURE);
	}

 	segmentArgs args;

	bool specifiedDetect = false;

 	/*parse the command line arguments */
	for ( int i = 1; i < argc; ){

 		std::string flag( argv[ i ] );

 		if ( flag == "-d" or flag == "--detect" ){
 			std::string strArg( argv[ i + 1 ] );
			args.detectFilename = strArg;
			i+=2;
			specifiedDetect = true;
		}
 		if ( flag == "-m" or flag == "--minLen" ){
			std::string strArg( argv[ i + 1 ] );
			args.minLength = std::stoi( strArg.c_str() );
			i+=2;
		}
		else if ( flag == "-t" or flag == "--threads" ){

			std::string strArg( argv[ i + 1 ] );
			args.threads = std::stoi( strArg.c_str() );
			i+=2;
		}
		else throw InvalidOption( flag );
	}

	if (args.minLength < 2){
 		std::cerr << "Exiting with error.  Minimum segment length must be at least 2." << std::endl;
 		std::cerr << help << std::endl;
		exit(EXIT_FAILURE);
	}

	if (args.minLength < 50){
 		std::cerr << "Warning: Low minimum segment lengths may lead to false positives." << std::endl;
 		std::cerr << "Use care in interpreting the results." << std::endl;
 		std::cerr << help << std::endl;
		exit(EXIT_FAILURE);
	}

	if (not (specifiedDetect)){
 		std::cerr << "Exiting with error.  Missing required arguments." << std::endl;
 		std::cerr << help << std::endl;
		exit(EXIT_FAILURE);
	}

	return args;
}


std::string writeSegmentHeader(segmentArgs &args, KMeansResult analougeIncorporation){

	auto t = std::time(nullptr);
	auto tm = *std::localtime(&t);
	std::ostringstream oss;
	oss << std::put_time(&tm, "%d/%m/%Y %H:%M:%S");
	auto str = oss.str();
	
	std::string compute = "CPU";

	std::string out;
	out += "#DetectFile " + args.detectFilename + "\n";
	out += "#Threads " + std::to_string(args.threads) + "\n";
	out += "#Compute " + compute + "\n";
	out += "#SystemStartTime " + str + "\n";
	out += "#Software " + std::string(getExePath()) + "\n";
	out += "#Version " + std::string(VERSION) + "\n";
	out += "#Commit " + std::string(getGitCommit()) + "\n";
	out += "#EstimatedRegionBrdU " + std::to_string(analougeIncorporation.centroid_1) + "\n";
	out += "#EstimatedRegionEdU " + std::to_string(analougeIncorporation.centroid_2) + "\n";

	return out;
}


void emptyBuffer(std::vector< DetectedRead > &buffer, segmentArgs args, std::ofstream &EdUFile, std::ofstream &BrdUFile, KMeansResult analogueIncorporation){

	#pragma omp parallel for schedule(dynamic) shared(args, analogueIncorporation) num_threads(args.threads)
	for ( auto b = buffer.begin(); b < buffer.end(); b++) {

		runDBSCAN(*b, analogueIncorporation, args.minLength);
		callSegmentation(*b, args.minLength);

		std::string BrdUOutput, EdUOutput;
		bool segmentToForks = false;

		std::pair<std::string,std::string> analogueOutputPair = writeAnalogueRegions(*b, segmentToForks);
		BrdUOutput = analogueOutputPair.first;
		EdUOutput = analogueOutputPair.second;

		#pragma omp critical
		{
			EdUFile << EdUOutput;
			BrdUFile << BrdUOutput;
		}
	}
	buffer.clear();
}


int segment_main( int argc, char** argv ){

	segmentArgs args = parseSegmentArguments( argc, argv );

	unsigned int maxBufferSize = 20*(args.threads);

	//get a read count
	int readCount = 0;
	std::string line;
	std::ifstream inFile( args.detectFilename );
	if ( not inFile.is_open() ) throw IOerror( args.detectFilename );
	while( std::getline( inFile, line ) ){

		if ( line.substr(0,1) == ">" ) readCount++;
	}	
	progressBar pb(readCount,true);
	inFile.close();
	
	//estimate analogue incorporation
	KMeansResult analogueIncorporation = estimateAnalogueIncorporation(args.detectFilename, readCount);

	//open all the files
	std::ofstream EdUFile, BrdUFile;

	BrdUFile.open("BrdU_DNAscent_segment.bed");
	BrdUFile << writeSegmentHeader(args, analogueIncorporation);
	if ( not BrdUFile.is_open() ) throw IOerror( "BrdU_DNAscent_segment.bed" );

	EdUFile.open("EdU_DNAscent_segment.bed");
	EdUFile << writeSegmentHeader(args, analogueIncorporation);
	if ( not EdUFile.is_open() ) throw IOerror( "EdU_DNAscent_segment.bed" );

	int failed = 0;

	std::vector< DetectedRead > readBuffer;
	int progress = 0;
	while( std::getline( inFile, line ) ){

		if ( line.substr(0,1) == "#"){
			continue;
		}
		else if ( line.substr(0,1) == ">" ){

			//check the read length on the back of the buffer
			if (readBuffer.size() > 0){

				if (readBuffer.back().positions.size() < 2000){
					failed++;
					readBuffer.pop_back();
				}
			}

			//empty the buffer if it's full
			if (readBuffer.size() >= maxBufferSize) emptyBuffer(readBuffer, args, EdUFile, BrdUFile, analogueIncorporation);

			progress++;
			pb.displayProgress( progress, failed, 0 );

			DetectedRead d;
			d.header = line;
			std::stringstream ssLine(line);
			std::string column;
			int cIndex = 0;
			while ( std::getline( ssLine, column, ' ' ) ){

				if ( cIndex == 0 ) d.readID = column;
				else if ( cIndex == 1 ) d.chromosome = column;
				else if ( cIndex == 2 ) d.mappingLower = std::stoi(column);
				else if ( cIndex == 3 ) d.mappingUpper = std::stoi(column);
				else if ( cIndex == 4 ) d.strand = column;
				else throw DetectParsing();
				cIndex++;
			}

			assert(d.mappingUpper > d.mappingLower);
			readBuffer.push_back(d);
		}
		else{

			std::string column;
			std::stringstream ssLine(line);
			int position = -1, cIndex = 0;
			AnalogueScore B, E;
			bool qualityOK = true;
			while ( std::getline( ssLine, column, '\t' ) ){

				if ( cIndex == 0 ){

					position = std::stoi(column);
				}
				else if ( cIndex == 1 ){

					E.set(std::stof(column));
				}
				else if ( cIndex == 2 ){

					B.set(std::stof(column));
				}
				else if ( cIndex == 3 ){
					assert(column == "*");
					qualityOK = false;
				}
				cIndex++;
			}

			readBuffer.back().alignmentQuality.push_back(qualityOK);
			readBuffer.back().positions.push_back(position);
			readBuffer.back().brduCalls.push_back(B.get());
			readBuffer.back().eduCalls.push_back(E.get());
		}
	}

	//empty the buffer at the end
	if (readBuffer.back().positions.size() < 2000){
		readBuffer.pop_back();
	}
	emptyBuffer(readBuffer, args, EdUFile, BrdUFile, analogueIncorporation);

	inFile.close();
	EdUFile.close();
	BrdUFile.close();

	std::cout << std::endl << "Done." << std::endl;

	return 0;
}

