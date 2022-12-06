//----------------------------------------------------------
// Copyright 2020 University of Cambridge
// Written by Michael A. Boemo (mb915@cam.ac.uk)
// This software is licensed under GPL-3.0.  You should have
// received a copy of the license with this software.  If
// not, please Email the author.
//----------------------------------------------------------

#ifndef SENSE_H
#define SENSE_H

#include <cassert>
#include <vector>
#include <iostream>
#include <fstream>
#include "error_handling.h"

int sense_main( int argc, char** argv );

struct ReadSegment{
	int leftmostCoord = 0;
	int leftmostIdx = 0;
	int rightmostCoord = 0;
	int rightmostIdx = 0;
	int partners = 0;
	double score = 0.0;
	std::vector<double> stress_signature;
};


struct KMeansResult{
	double centroid_1;
	double centroid_1_lowerBound;
	double centroid_1_stdv;
	double centroid_2;
	double centroid_2_lowerBound;
	double centroid_2_stdv;
};

struct forkSenseArgs {

	std::string detectFilename;
	std::string outputFilename;
	std::string analogueOrder;
	bool deprecatedMarks = false;
	bool makeSignatures = false;
	unsigned int threads = 1;
};

class AnalogueScore{

	private:
		double _score = 0.0;
		bool _isSet = false;
	public:
		void set(double s){
			_score = s;
			_isSet = true;
		}
		double get(void){
			assert(_isSet);
			return _score;
		}
};

class DetectedRead{

	public:
		std::vector< int > positions;
		std::vector< bool > alignmentQuality;
		std::vector< double > brduCalls, eduCalls;
		std::map< int, double > stallScore;
		std::string readID, chromosome, strand, header;
		int mappingLower, mappingUpper;
		std::vector< int > EdU_segment_label, BrdU_segment_label, thymidine_segment_label;
		std::vector<ReadSegment> EdU_segment, BrdU_segment;
		std::vector<ReadSegment> origins, terminations, leftForks, rightForks;
		std::vector<double> tensorInput;
		int64_t inputSize;
};


std::string writeForkSenseHeader(forkSenseArgs &, KMeansResult );
std::string writeBedHeader( forkSenseArgs & );

class fs_fileManager{

	protected:
	
		forkSenseArgs inputArgs;
		std::ofstream outFile, originFile, termFile, leftForkFile, rightForkFile, leftSignaturesFile, rightSignaturesFile, EdUFile, BrdUFile;

	public:
	
		fs_fileManager( forkSenseArgs &args , KMeansResult analogueIncorporation){
		
			inputArgs = args;
		
			//main output file
		 	outFile.open( args.outputFilename );
			if ( not outFile.is_open() ) throw IOerror( args.outputFilename );
			std::string outHeader = writeForkSenseHeader(args, analogueIncorporation);
			outFile << outHeader;
			
			//bed files
			termFile.open("terminations_DNAscent_forkSense.bed");
			termFile << writeBedHeader(args);
			if ( not termFile.is_open() ) throw IOerror( "terminations_DNAscent_forkSense.bed" );

			originFile.open("origins_DNAscent_forkSense.bed");
			originFile << writeBedHeader(args);
			if ( not originFile.is_open() ) throw IOerror( "origins_DNAscent_forkSense.bed" );

			leftForkFile.open("leftForks_DNAscent_forkSense.bed");
			leftForkFile << writeBedHeader(args);
			if ( not leftForkFile.is_open() ) throw IOerror( "leftForks_DNAscent_forkSense.bed" );

			rightForkFile.open("rightForks_DNAscent_forkSense.bed");
			rightForkFile << writeBedHeader(args);
			if ( not rightForkFile.is_open() ) throw IOerror( "rightForks_DNAscent_forkSense.bed" );

			BrdUFile.open("BrdU_DNAscent_forkSense.bed");
			BrdUFile << writeBedHeader(args);
			if ( not BrdUFile.is_open() ) throw IOerror( "BrdU_DNAscent_forkSense.bed" );

			EdUFile.open("EdU_DNAscent_forkSense.bed");
			EdUFile << writeBedHeader(args);
			if ( not EdUFile.is_open() ) throw IOerror( "EdU_DNAscent_forkSense.bed" );

			if (args.makeSignatures){

				leftSignaturesFile.open("leftForks_DNAscent_forkSense_stressSignatures.bed");
				leftSignaturesFile << writeBedHeader(args);
				if ( not leftSignaturesFile.is_open() ) throw IOerror( "leftForks_DNAscent_forkSense_stressSignatures.bed" );

				rightSignaturesFile.open("rightForks_DNAscent_forkSense_stressSignatures.bed");
				rightSignaturesFile << writeBedHeader(args);
				if ( not rightSignaturesFile.is_open() ) throw IOerror( "rightForks_DNAscent_forkSense_stressSignatures.bed" );
			}
		}
		void writeOutput(std::string &readOutput,
				std::string &termOutput,
				std::string &originOutput,		
				std::string &leftForkOutput,		
				std::string &rightForkOutput,	
				std::string &leftSignaturesOutput,		
				std::string &rightSignaturesOutput,		
				std::string &BrdUOutput,		
				std::string &EdUOutput	){
		
			outFile << readOutput;
			termFile << termOutput;
			originFile << originOutput;
			leftForkFile << leftForkOutput;
			rightForkFile << rightForkOutput;
			BrdUFile << BrdUOutput;
			EdUFile << EdUOutput;

			if (inputArgs.makeSignatures and leftSignaturesOutput.size() > 0) leftSignaturesFile << leftSignaturesOutput;
			if (inputArgs.makeSignatures and rightSignaturesOutput.size() > 0) rightSignaturesFile << rightSignaturesOutput;

		}
		void closeAll(){
			outFile.close();
			termFile.close();
			originFile.close();
			BrdUFile.close();
			EdUFile.close();
			leftForkFile.close();
			rightForkFile.close();

			if (inputArgs.makeSignatures){
				leftSignaturesFile.close();
				rightSignaturesFile.close();			
			}
		}
};

KMeansResult twoMeans_fs( std::vector< double > & );
std::pair<int, int> segmentationTrim(std::vector< int > &, std::vector< double > &, std::vector< double > &, int , int );
KMeansResult estimateAnalogueIncorporation(std::string , int );
void runDBSCAN(DetectedRead &, KMeansResult, int, double);
void callSegmentation(DetectedRead &, int);
std::pair<std::string,std::string> writeAnalogueRegions(DetectedRead &, bool );


#endif

