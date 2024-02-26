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
#include "tensor.h"
#include <cmath>
#include <memory>
#include <math.h>
#include <algorithm>
#include <limits>
#include <stdlib.h>

static const char *help=
"forkSense: DNAscent AI executable that calls replication origins, termination sites, and fork movement.\n"
"To run DNAscent forkSense, do:\n"
"   DNAscent forkSense -d /path/to/BrdUCalls.detect -o /path/to/output.forkSense --order EdU,BrdU\n"
"Required arguments are:\n"
"  -d,--detect               path to output file from DNAscent detect,\n"
"  -o,--output               path to output file for forkSense,\n"
"     --order                order in which the analogues were pulsed (EdU,BrdU or BrdU,EdU).\n"
"Optional arguments are:\n"
"  -t,--threads              number of threads (default: 1 thread),\n"
"  --markAnalogues           writes analogue incorporation locations to a bed file (default: off),\n"
"  --markOrigins             writes replication origin locations to a bed file (default: off),\n"
"  --markTerminations        writes replication termination locations to a bed file (default: off),\n"
"  --markForks               writes replication fork locations to a bed file (default: off),\n"
"  --makeSignatures          writes replication stress signatures to a bed files (default: off).\n"
"DNAscent is under active development by the Boemo Group, Department of Pathology, University of Cambridge (https://www.boemogroup.org/).\n"
"Please submit bug reports to GitHub Issues (https://github.com/MBoemo/DNAscent/issues).";


forkSenseArgs parseSenseArguments( int argc, char** argv ){

 	if( argc < 2 ){
 		std::cout << "Exiting with error.  Insufficient arguments passed to DNAscent forkSense." << std::endl << help << std::endl;
		exit(EXIT_FAILURE);
	}
 	if ( std::string( argv[ 1 ] ) == "-h" or std::string( argv[ 1 ] ) == "--help" ){
 		std::cout << help << std::endl;
		exit(EXIT_SUCCESS);
	}
	else if( argc < 4 ){ //PLP&SY: check with Mike
 		std::cout << "Exiting with error.  Insufficient arguments passed to DNAscent forkSense." << std::endl;
		exit(EXIT_FAILURE);
	}

 	forkSenseArgs args;

	bool specifiedDetect = false;
	bool specifiedOutput = false;
	bool specifiedOrder = false;

 	/*parse the command line arguments */
	for ( int i = 1; i < argc; ){

 		std::string flag( argv[ i ] );

 		if ( flag == "-d" or flag == "--detect" ){
 			std::string strArg( argv[ i + 1 ] );
			args.detectFilename = strArg;
			i+=2;
			specifiedDetect = true;
		}
		else if ( flag == "-o" or flag == "--output" ){
 			std::string strArg( argv[ i + 1 ] );
			args.outputFilename = strArg;
			i+=2;
			specifiedOutput = true;
		}
		else if (flag == "--order" ){
 			std::string strArg( argv[ i + 1 ] );
			args.analogueOrder = strArg;
			i+=2;
			specifiedOrder = true;
		}
		else if ( flag == "-t" or flag == "--threads" ){

			std::string strArg( argv[ i + 1 ] );
			args.threads = std::stoi( strArg.c_str() );
			i+=2;
		}
		else if ( flag == "--markOrigins" ){

			args.markOrigins = true;
			i+=1;
		}
		else if ( flag == "--markTerminations" ){

			args.markTerms = true;
			i+=1;
		}
		else if ( flag == "--markForks" ){

			args.markForks = true;
			i+=1;
		}
		else if ( flag == "--markAnalogues" ){

			args.markAnalogues = true;
			i+=1;
		}
		else if ( flag == "--makeSignatures" ){

			args.makeSignatures = true;
			i+=1;
		}
		else throw InvalidOption( flag );
	}
	if (args.outputFilename == args.detectFilename) throw OverwriteFailure();
	if (not (specifiedDetect and specifiedOutput and specifiedOrder)){
 		std::cout << "Exiting with error.  Missing required arguments." << std::endl;
 		std::cout << help << std::endl;
		exit(EXIT_FAILURE);
	}

	if (args.analogueOrder != "EdU,BrdU" and args.analogueOrder != "BrdU,EdU"){
 		std::cout << "Exiting with error.  Analogue order should be EdU,BrdU or BrdU,EdU." << std::endl;
 		std::cout << help << std::endl;
		exit(EXIT_FAILURE);
	}

	return args;
}


std::string writeForkSenseHeader(forkSenseArgs &args, KMeansResult analougeIncorporation){

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


std::string writeBedHeader( forkSenseArgs &args){

	auto t = std::time(nullptr);
	auto tm = *std::localtime(&t);
	std::ostringstream oss;
	oss << std::put_time(&tm, "%d/%m/%Y %H:%M:%S");
	auto str = oss.str();
	
	std::string compute = "CPU";

	std::string out;
	out += "#DetectFile " + args.detectFilename + "\n";
	out += "#ForkSenseFile " + args.outputFilename + "\n";
	out += "#AnalogueOrder " + args.analogueOrder + "\n";
	out += "#Threads " + std::to_string(args.threads) + "\n";
	out += "#Compute " + compute + "\n";
	out += "#SystemStartTime " + str + "\n";
	out += "#Software " + std::string(getExePath()) + "\n";
	out += "#Version " + std::string(VERSION) + "\n";
	out += "#Commit " + std::string(getGitCommit()) + "\n";

	return out;
}


std::vector<ReadSegment> stitchSegmentation(std::vector<ReadSegment> &primarySegments, std::vector<ReadSegment> &secondarySegments){

	std::map<int, int> connectivity;
	std::vector<ReadSegment> stitchedSegments;
	
	int segmentStitch = 3000;
	
	for (size_t i = 0; i < primarySegments.size(); i++){
	
		for (size_t j = i+1; j < primarySegments.size(); j++){
		
			//segments shouldn't overlap
			assert(primarySegments[j].leftmostCoord >= primarySegments[i].rightmostCoord or primarySegments[i].leftmostCoord >= primarySegments[j].rightmostCoord);
		
			if (primarySegments[j].leftmostCoord - primarySegments[i].rightmostCoord < segmentStitch){
			
				//make sure there's not a short segment between
				bool interveningSegment = false;
				for (size_t k = 0; k < secondarySegments.size(); k++){
				
					//segments shouldn't overlap
					assert(primarySegments[i].leftmostCoord >= secondarySegments[k].rightmostCoord or secondarySegments[k].leftmostCoord >= primarySegments[i].rightmostCoord);
					
					if (primarySegments[i].rightmostCoord <= secondarySegments[k].leftmostCoord and secondarySegments[k].rightmostCoord <= primarySegments[j].leftmostCoord){
						interveningSegment = true;
						break;
					}
				}
				
				if (not interveningSegment){
					connectivity[i] = j;
					break;
				}
			}
		}
	}
	
	//stitch 
	std::vector<int> ignoreIndices;
	for (size_t i = 0; i < primarySegments.size(); i++){
	
		if (std::find(ignoreIndices.begin(), ignoreIndices.end(), i) != ignoreIndices.end()) continue;
	
		int targetLeft = i;
		int startCoord = primarySegments[targetLeft].leftmostCoord;
		int startIdx = primarySegments[targetLeft].leftmostIdx;
		int endCoord = primarySegments[targetLeft].rightmostCoord;
		int endIdx = primarySegments[targetLeft].rightmostIdx;
		
		while (connectivity.count(targetLeft) > 0){
		
			int idxToMerge = connectivity[targetLeft];
		
			assert(endCoord < primarySegments[idxToMerge].rightmostCoord);
			assert(endIdx < primarySegments[idxToMerge].rightmostIdx);
		
			endCoord = primarySegments[idxToMerge].rightmostCoord;
			endIdx = primarySegments[idxToMerge].rightmostIdx;
			ignoreIndices.push_back(idxToMerge);
			targetLeft = idxToMerge;
		}
	
		struct ReadSegment s = {startCoord, startIdx, endCoord, endIdx};
		stitchedSegments.push_back(s);
	}
	return stitchedSegments;
}


void callSegmentation(DetectedRead &r){

	int minLength = 1000;

	bool inSegment = false;
	int startCoord = -1, endCoord = -1;
	int startIdx = -1, endIdx = -1;
	
	std::vector<ReadSegment> EdU_segments;
	std::vector<ReadSegment> BrdU_segments;

	for (size_t i = 0; i < r.positions.size(); i++){

		if (r.EdU_segment_label[i] == 1 and not inSegment){ //initialise the site

			startCoord = r.positions[i];
			startIdx = i;
			inSegment = true;
		}
		else if (inSegment and (r.thymidine_segment_label[i] == 1 or r.BrdU_segment_label[i] == 1)){//close if we've confidently moved to something else
		
			endCoord = r.positions[i];
			endIdx = i;

			assert(startCoord != -1 and endCoord != -1);

			if ( abs(endCoord - startCoord) >= minLength ){

				//std::pair<int, int> trim = segmentationTrim(r.positions, r.eduCalls, r.brduCalls, startIdx, endIdx);
				//startIdx += trim.first;
				//endIdx -= trim.second;
				startCoord = r.positions[startIdx];
				endCoord = r.positions[endIdx];
				
				assert(startCoord < endCoord);
				
				struct ReadSegment s = {startCoord, startIdx, endCoord, endIdx};
				EdU_segments.push_back(s);
			}

			inSegment = false;
			startCoord = -1;
			endCoord = -1;
			startIdx = -1;
			endIdx = -1;
		}
	}
	//if we got to the end of the read without closing
	if (inSegment){

		assert(startCoord != -1);
		if (endCoord == -1){
			endCoord = r.positions.back();
			endIdx = r.positions.size() - 1;
		}

		if ( abs(endCoord - startCoord) >= minLength ){

			//std::pair<int, int> trim = segmentationTrim(r.positions, r.eduCalls, r.brduCalls, startIdx, endIdx);
			//startIdx += trim.first;
			//endIdx -= trim.second;
			startCoord = r.positions[startIdx];
			endCoord = r.positions[endIdx];
			
			assert(startCoord < endCoord);
			
			struct ReadSegment s = {startCoord, startIdx, endCoord, endIdx};
			EdU_segments.push_back(s);
		}
	}

	startCoord = -1;
	endCoord = -1;
	startIdx = -1;
	endIdx = -1;
	inSegment = false;

	for (size_t i = 0; i < r.positions.size(); i++){

		if (r.BrdU_segment_label[i] == 1 and not inSegment){ //initialise the site

			startCoord = r.positions[i];
			startIdx = i;
			inSegment = true;
		}
		else if (inSegment and (r.thymidine_segment_label[i] == 1 or r.EdU_segment_label[i] == 1)){//close if we've confidently moved to something else
		
			endCoord = r.positions[i];
			endIdx = i;

			assert(startCoord != -1 and endCoord != -1);

			if ( abs(endCoord - startCoord) >= minLength ){
			
				//std::pair<int, int> trim = segmentationTrim(r.positions, r.brduCalls, r.eduCalls, startIdx, endIdx);
				//startIdx += trim.first;
				//endIdx -= trim.second;
				startCoord = r.positions[startIdx];
				endCoord = r.positions[endIdx];
				
				assert(startCoord < endCoord);
			
				struct ReadSegment s = {startCoord, startIdx, endCoord, endIdx};
				BrdU_segments.push_back(s);
			}

			inSegment = false;
			startCoord = -1;
			endCoord = -1;
			startIdx = -1;
			endIdx = -1;
		}
	}
	//if we got to the end of the read without closing
	if (inSegment){

		assert(startCoord != -1);
		if (endCoord == -1){
			endCoord = r.positions.back();
			endIdx = r.positions.size() - 1;
		}

		if ( abs(endCoord - startCoord) >= minLength ){
		
			//std::pair<int, int> trim = segmentationTrim(r.positions, r.brduCalls, r.eduCalls, startIdx, endIdx);
			//startIdx += trim.first;
			//endIdx -= trim.second;
			startCoord = r.positions[startIdx];
			endCoord = r.positions[endIdx];
			
			assert(startCoord < endCoord);
		
			struct ReadSegment s = {startCoord, startIdx, endCoord, endIdx};
			BrdU_segments.push_back(s);
		}
	}
	
	r.BrdU_segment = stitchSegmentation(BrdU_segments, EdU_segments);
	r.EdU_segment = stitchSegmentation(EdU_segments, BrdU_segments);
}


std::string callOrigins(DetectedRead &r, forkSenseArgs args){

	std::string outBed;

	//match up regions
	for ( size_t li = 0; li < r.leftForks.size(); li++ ){

		//find the closest right fork region
		int minDist = std::numeric_limits<int>::max();
		int bestMatch = -1;
		for ( size_t ri = 0; ri < r.rightForks.size(); ri++ ){

			if (r.rightForks[ri].rightmostCoord < r.leftForks[li].rightmostCoord) continue;

			int dist = r.rightForks[ri].rightmostCoord - r.leftForks[li].leftmostCoord;
			assert(dist >= 0);
			if (dist < minDist){
				minDist = dist;
				bestMatch = ri;
			}
		}

		//make sure no other left forks are closer
		bool failed = false;
		if (bestMatch != -1){

			for (size_t l2 = 0; l2 < r.leftForks.size(); l2++){

				if (l2 == li) continue;

				if (r.rightForks[bestMatch].rightmostCoord < r.leftForks[l2].rightmostCoord) continue;

				int dist = r.rightForks[bestMatch].rightmostCoord - r.leftForks[l2].leftmostCoord;
				assert(dist >= 0);

				if (dist < minDist){

					failed = true;
					break;
				}
			}
		}
		if (failed) continue;
		else if (bestMatch != -1){

			int orilb = std::min(r.leftForks[li].rightmostCoord, r.rightForks[bestMatch].leftmostCoord);
			int oriub = std::max(r.leftForks[li].rightmostCoord, r.rightForks[bestMatch].leftmostCoord);
			
			int orilb_idx = std::min(r.leftForks[li].rightmostIdx, r.rightForks[bestMatch].leftmostIdx);
			int oriub_idx = std::max(r.leftForks[li].rightmostIdx, r.rightForks[bestMatch].leftmostIdx);
			
			struct ReadSegment s = {orilb, orilb_idx, oriub, oriub_idx};
			r.origins.push_back(s);

			outBed += r.chromosome + " " 
			       + std::to_string(orilb) + " " 
			       + std::to_string(oriub) + " "
			       + r.readID.substr(1) + " "
			       + std::to_string( r.mappingLower ) + " "
			       + std::to_string( r.mappingUpper ) + " "
			       + r.strand + "\n"; 
		}
	}

	return outBed;
}


std::string callTerminations(DetectedRead &r, forkSenseArgs args){

	std::string outBed;

	//match up regions
	for ( size_t li = 0; li < r.leftForks.size(); li++ ){

		//find the closest right fork region
		int minDist = std::numeric_limits<int>::max();
		int bestMatch = -1;
		for ( size_t ri = 0; ri < r.rightForks.size(); ri++ ){

			if (r.leftForks[li].rightmostCoord < r.rightForks[ri].rightmostCoord ) continue;

			int dist = r.leftForks[li].rightmostCoord - r.rightForks[ri].leftmostCoord;
			assert(dist >= 0);

			if (dist < minDist){
				minDist = dist;
				bestMatch = ri;

			}
		}

		//make sure no other left forks are closer
		bool failed = false;
		if (bestMatch != -1){

			for (size_t l2 = 0; l2 < r.leftForks.size(); l2++){

				if (li == l2) continue;

				if (r.leftForks[l2].rightmostCoord < r.rightForks[bestMatch].rightmostCoord ) continue;

				int dist = r.leftForks[l2].rightmostCoord - r.rightForks[bestMatch].leftmostCoord;
				assert(dist >= 0);

				if (dist < minDist){

					failed = true;
					break;
				}
			}
		}
		if (failed) continue;
		else if (bestMatch != -1){

			int termlb = std::min(r.leftForks[li].leftmostCoord,r.rightForks[bestMatch].rightmostCoord);
			int termub = std::max(r.leftForks[li].leftmostCoord,r.rightForks[bestMatch].rightmostCoord);
			
			int termlb_idx = std::min(r.leftForks[li].leftmostIdx,r.rightForks[bestMatch].rightmostIdx);
			int termub_idx = std::max(r.leftForks[li].leftmostIdx,r.rightForks[bestMatch].rightmostIdx);
			
			struct ReadSegment s = {termlb, termlb_idx, termub, termub_idx};

			r.terminations.push_back(s);
			outBed += r.chromosome + " " 
			       + std::to_string(termlb) + " " 
			       + std::to_string(termub) + " " 
			       + r.readID.substr(1) + " "
			       + std::to_string( r.mappingLower ) + " "
			       + std::to_string( r.mappingUpper ) + " "
			       + r.strand + "\n"; 
		}
	}

	return outBed;
}


std::pair<std::string,std::string> writeAnalogueRegions(DetectedRead &r, bool segmentToForks){

	std::string outBedEdU, outBedBrdU;

	for (auto i = r.BrdU_segment.begin(); i < r.BrdU_segment.end(); i++){
	
		if (segmentToForks and (i -> partners == 0)) continue;

		outBedBrdU += r.chromosome + " " 
		            + std::to_string(i -> leftmostCoord) + " " 
		            + std::to_string(i -> rightmostCoord) + " " 
		            + r.readID.substr(1) + " "
		            + std::to_string( r.mappingLower ) + " "
		            + std::to_string( r.mappingUpper ) + " "
		            + r.strand + "\n"; 
	}
	for (auto i = r.EdU_segment.begin(); i < r.EdU_segment.end(); i++){
	
		if (segmentToForks and (i -> partners == 0)) continue;

		outBedEdU += r.chromosome + " " 
		            + std::to_string(i -> leftmostCoord) + " " 
		            + std::to_string(i -> rightmostCoord) + " " 
		            + r.readID.substr(1) + " "
		            + std::to_string( r.mappingLower ) + " "
		            + std::to_string( r.mappingUpper ) + " "
		            + r.strand + "\n"; 
	}

	return std::make_pair(outBedBrdU,outBedEdU);
}


void callForks(DetectedRead &r, std::string analogueOrder){

	//maximum distance in bp between EdU and BrdU segments that allow them to be matched
	int maxGap = 5000;
	
	std::vector<ReadSegment> *analogue1_segments;
	std::vector<ReadSegment> *analogue2_segments;
	
	if (analogueOrder == "EdU,BrdU"){

		analogue1_segments = &(r.EdU_segment);
		analogue2_segments = &(r.BrdU_segment);
	}
	else{

		analogue2_segments = &(r.EdU_segment);
		analogue1_segments = &(r.BrdU_segment);
	}

	std::vector<std::pair<size_t,size_t>> proto_rightForkPairs, proto_leftForkPairs;

	//match up segments - right fork
	for ( size_t li = 0; li < (*analogue1_segments).size(); li++ ){

		//find the closest analogue 2 patch
		int minDist = std::numeric_limits<int>::max();
		int bestMatch = -1;
		for ( size_t ri = 0; ri < (*analogue2_segments).size(); ri++ ){

			if ( (*analogue2_segments)[ri].leftmostCoord < (*analogue1_segments)[li].rightmostCoord ) continue;

			int dist = (*analogue2_segments)[ri].leftmostCoord - (*analogue1_segments)[li].rightmostCoord;
			assert(dist >= 0);
			
			if (dist < minDist){
				minDist = dist;
				bestMatch = ri;
			}
		}

		//make sure no other analogue 1 patches are closer
		bool failed = false;
		if (bestMatch != -1){

			for (size_t l2 = 0; l2 < (*analogue1_segments).size(); l2++){

				if (li == l2) continue;

				if ( (*analogue2_segments)[bestMatch].leftmostCoord < (*analogue1_segments)[l2].rightmostCoord ) continue;

				int dist = (*analogue2_segments)[bestMatch].leftmostCoord - (*analogue1_segments)[l2].rightmostCoord;
				assert(dist >= 0);
				
				if (dist < minDist){

					failed = true;
					break;
				}
			}
		}

		if (failed) continue;
		else if (bestMatch != -1 and minDist < maxGap){

			assert( (*analogue1_segments)[li].leftmostCoord < (*analogue2_segments)[bestMatch].rightmostCoord );
			
			(*analogue1_segments)[li].partners++;
			(*analogue2_segments)[bestMatch].partners++;
			
			proto_rightForkPairs.push_back(std::make_pair(li,bestMatch));
		}
	}


	//match up segments - left fork
	for ( size_t li = 0; li < (*analogue1_segments).size(); li++ ){

		//find the closest analogue 2 patch
		int minDist = std::numeric_limits<int>::max();
		int bestMatch = -1;
		for ( size_t ri = 0; ri < (*analogue2_segments).size(); ri++ ){

			if ( (*analogue1_segments)[li].leftmostCoord < (*analogue2_segments)[ri].rightmostCoord ) continue;
			
			int dist = (*analogue1_segments)[li].leftmostCoord - (*analogue2_segments)[ri].rightmostCoord;
			assert(dist >= 0);

			if (dist < minDist){
				minDist = dist;
				bestMatch = ri;
			}
		}

		//make sure no other analogue 2 patches are closer
		bool failed = false;
		if (bestMatch != -1){

			for (size_t l2 = 0; l2 < (*analogue1_segments).size(); l2++){

				if (li == l2) continue;

				if ( (*analogue1_segments)[l2].leftmostCoord <  (*analogue2_segments)[bestMatch].rightmostCoord ) continue;

				int dist = (*analogue1_segments)[l2].leftmostCoord -  (*analogue2_segments)[bestMatch].rightmostCoord;
				assert(dist >= 0);
				
				if (dist < minDist){
					failed = true;
					break;
				}
			}
		}
		if (failed) continue;
		else if (bestMatch != -1 and minDist < maxGap){

			assert( (*analogue2_segments)[bestMatch].leftmostCoord < (*analogue1_segments)[li].rightmostCoord );
			
			(*analogue2_segments)[bestMatch].partners++;
			(*analogue1_segments)[li].partners++;
			
			proto_leftForkPairs.push_back(std::make_pair(bestMatch,li));
		}
	}
	
	//make fork bounds
	std::string outBedLeft, outBedRight;
	for (auto p = proto_rightForkPairs.begin(); p < proto_rightForkPairs.end(); p++){
	
		int tipPartners = 0;
	
		//left coordinate and index
		assert( ((*analogue1_segments)[p -> first].partners == 1) or ((*analogue1_segments)[p -> first].partners == 2));
		int lc = (*analogue1_segments)[p -> first].leftmostCoord;
		int li = (*analogue1_segments)[p -> first].leftmostIdx;
		if ( (*analogue1_segments)[p -> first].partners == 2){
		
			lc = ( (*analogue1_segments)[p -> first].leftmostCoord + (*analogue1_segments)[p -> first].rightmostCoord ) / 2;
			li = ( (*analogue1_segments)[p -> first].leftmostIdx + (*analogue1_segments)[p -> first].rightmostIdx ) / 2;
		}

		//right coordinate and index		
		int rc = (*analogue2_segments)[p -> second].rightmostCoord;
		int ri = (*analogue2_segments)[p -> second].rightmostIdx;
		if ( (*analogue2_segments)[p -> second].partners == 2){
		
			rc = ( (*analogue2_segments)[p -> second].rightmostCoord + (*analogue2_segments)[p -> second].leftmostCoord ) / 2;
			ri = ( (*analogue2_segments)[p -> second].rightmostIdx + (*analogue2_segments)[p -> second].leftmostIdx ) / 2;
			tipPartners++;
		}
		
		double An1Len = (*analogue1_segments)[p -> first].rightmostCoord - lc;
		double An2Len = rc - (*analogue2_segments)[p -> second].leftmostCoord;
		
		int BrdU_in_EdU = 0;
		int EdU_in_EdU = 0;
		int attempts_EdU = 0;
		for (int j = li; j < (*analogue1_segments)[p -> first].rightmostIdx; j++){
		
			if (r.brduCalls[j] > 0.5) BrdU_in_EdU++;
			if (r.eduCalls[j] > 0.5) EdU_in_EdU++;
			attempts_EdU++;
		}

		int BrdU_in_BrdU = 0;
		int EdU_in_BrdU = 0;
		int attempts_BrdU = 0;	
		for (int j = (*analogue2_segments)[p -> second].leftmostIdx; j < ri; j++){
		
			if (r.brduCalls[j] > 0.5) BrdU_in_BrdU++;
			if (r.eduCalls[j] > 0.5) EdU_in_BrdU++;
			attempts_BrdU++;
		}
		
		struct ReadSegment s = {lc, li, rc, ri};
		s.partners = tipPartners;
		
		double overallLength = rc-lc;
		double sig4 = BrdU_in_EdU/(double) attempts_EdU;
		double sig5 = EdU_in_EdU/(double) attempts_EdU;
		double sig6 = EdU_in_BrdU/(double) attempts_BrdU;
		double sig7 = BrdU_in_BrdU/(double) attempts_BrdU;
		s.stress_signature = {overallLength,
					An1Len,
					An2Len, 
					sig4, 
					sig5, 
					sig6, 
					sig7};
		r.rightForks.push_back(s);
	}
	for (auto p = proto_leftForkPairs.begin(); p < proto_leftForkPairs.end(); p++){
	
		int tipPartners = 0;

		//left coordinate and index
		assert( ((*analogue2_segments)[p -> first].partners == 1) or ((*analogue2_segments)[p -> first].partners == 2));
		int lc = (*analogue2_segments)[p -> first].leftmostCoord;
		int li = (*analogue2_segments)[p -> first].leftmostIdx;
		if ( (*analogue2_segments)[p -> first].partners == 2){
		
			lc = ( (*analogue2_segments)[p -> first].leftmostCoord + (*analogue2_segments)[p -> first].rightmostCoord ) / 2;
			li = ( (*analogue2_segments)[p -> first].leftmostIdx + (*analogue2_segments)[p -> first].rightmostIdx ) / 2;
			tipPartners++;
		}
		
		//right coordinate and index
		int rc = (*analogue1_segments)[p -> second].rightmostCoord;
		int ri = (*analogue1_segments)[p -> second].rightmostIdx;
		if ( (*analogue1_segments)[p -> second].partners == 2){
		
			rc = ( (*analogue1_segments)[p -> second].rightmostCoord + (*analogue1_segments)[p -> second].leftmostCoord ) / 2;
			ri = ( (*analogue1_segments)[p -> second].rightmostIdx + (*analogue1_segments)[p -> second].leftmostIdx ) / 2;
		}
		
		double An2Len = (*analogue2_segments)[p -> first].rightmostCoord - lc;
		double An1Len = rc - (*analogue1_segments)[p -> second].leftmostCoord;

		int BrdU_in_EdU = 0;
		int EdU_in_EdU = 0;
		int attempts_EdU = 0;
		
		for (int j = (*analogue1_segments)[p -> second].leftmostIdx; j < ri; j++){
		
			if (r.brduCalls[j] > 0.5) BrdU_in_EdU++;
			if (r.eduCalls[j] > 0.5) EdU_in_EdU++;
			attempts_EdU++;
		}
		
		int BrdU_in_BrdU = 0;
		int EdU_in_BrdU = 0;
		int attempts_BrdU = 0;
		
		for (int j = li; j < (*analogue2_segments)[p -> first].rightmostIdx; j++){
		
			if (r.brduCalls[j] > 0.5) BrdU_in_BrdU++;
			if (r.eduCalls[j] > 0.5) EdU_in_BrdU++;
			attempts_BrdU++;
		}
		
		struct ReadSegment s = {lc, li, rc, ri};
		
		double overallLength = rc-lc;
		double sig4 = BrdU_in_EdU/(double) attempts_EdU;
		double sig5 = EdU_in_EdU/(double) attempts_EdU;
		double sig6 = EdU_in_BrdU/(double) attempts_BrdU;
		double sig7 = BrdU_in_BrdU/(double) attempts_BrdU;
		s.stress_signature = {overallLength,
					An1Len,
					An2Len, 
					sig4, 
					sig5, 
					sig6, 
					sig7};
					
		s.partners = tipPartners;
		r.leftForks.push_back(s);
	}
}


std::pair< std::vector<int>, int > findNeighbours_mod( std::vector<int> &positions, std::vector< double > &calls, std::vector< double > &altCalls, int index, int epsilon ){

	std::vector< int > neighbourIdx;
	int positiveCalls = 0;
	int positiveAltCalls = 0;
	int ev = positions[index];
	
	int windowStart = index-epsilon;
	int windowEnd = index+epsilon;
	int numPositions = positions.size();
	int startIdx = std::max(windowStart, 0);
	int endIdx = std::min(windowEnd, numPositions-1);
	
	for ( int i = startIdx; i <= endIdx; i++ ){
		int runningPos = positions[i];
		int gap = std::abs(ev - runningPos);

		if (gap <= epsilon){
		
			neighbourIdx.push_back(i);

			if (calls[i] > 0.5){
				positiveCalls++;
			}
			if (altCalls[i] > 0.5){
				positiveAltCalls++;
			}
		}
	}
	
	int delta = positiveCalls - positiveAltCalls;
	int netPositiveCalls = std::max(0,delta);

	
	return std::make_pair(neighbourIdx, netPositiveCalls);
}

std::map<int,int> DBSCAN_mod( std::vector< int > &positions, std::vector< double > &calls, std::vector< double > &altCalls, int epsilon, double minDensity ){

	//initialise labels
	std::map< int, int > index2label;
	for ( size_t i = 0; i < positions.size(); i++ ) index2label[i] = -2;

	for ( size_t i = 0; i < positions.size(); i++ ){
	
		std::pair<std::vector<int>, int> neighbourPair = findNeighbours_mod( positions, calls, altCalls, i, epsilon );
		std::vector<int> neighbourIndices = neighbourPair.first;
		int neighbourCalls = neighbourPair.second;
		int minPoints = neighbourIndices.size() * minDensity;
		if (neighbourCalls < minPoints) {

			index2label[i] = -1; //label as noise
		}
		else{
		
			index2label[i] = 1;
		}
	}
	return index2label;
}


void runDBSCAN(DetectedRead &r, KMeansResult analougeIncorporation, std::string analogueOrder){
	
	int epsilon = 500;
	
	double minBrdUDensity, minEdUDensity;
	
	if (analogueOrder == "EdU,BrdU"){

		minBrdUDensity = std::max(0.1, analougeIncorporation.centroid_1_lowerBound);
		minEdUDensity = std::max(0.1, analougeIncorporation.centroid_2_lowerBound);
	}
	else{

		minBrdUDensity = std::max(0.1, analougeIncorporation.centroid_1_lowerBound);
		minEdUDensity = std::max(0.1, analougeIncorporation.centroid_2_lowerBound);
	}
	
	
	std::map<int,int> eduLabels = DBSCAN_mod( r.positions, r.eduCalls, r.brduCalls, epsilon, minEdUDensity );
	std::map<int,int> brduLabels = DBSCAN_mod( r.positions, r.brduCalls, r.eduCalls, epsilon, minBrdUDensity );
	

	for (size_t i = 0; i < r.positions.size(); i++){
	
		int eduLabel = 0;
		int brduLabel = 0;
		int thymLabel = 0;
		
		if (eduLabels[ i ] >= 0 and brduLabels[ i ] < 0){
			eduLabel = 1;
			brduLabel = 0;
			thymLabel = 0;
		}	
		else if (brduLabels[ i ] >= 0 and eduLabels[ i ] < 0){
			eduLabel = 0;
			brduLabel = 1;
			thymLabel = 0;
		}
		else if (brduLabels[ i ] < 0 and eduLabels[ i ] < 0){
			eduLabel = 0;
			brduLabel = 0;
			thymLabel = 1;
		}

		r.EdU_segment_label.push_back(eduLabel);
		r.BrdU_segment_label.push_back(brduLabel);
		r.thymidine_segment_label.push_back(thymLabel);
	}
}


std::pair<int, int> segmentationTrim(std::vector< int > &positions, std::vector< double > &calls, std::vector< double > &altCalls, int startIdx, int endIdx){

	int epsilon = 500; //500 bp window
	
	if (positions[endIdx] - positions[startIdx] < 3*epsilon){
	
		return std::make_pair(0,0);
	}
	
	std::vector< int > segmentPositions(positions.begin() + startIdx, positions.begin() + endIdx + 1);
	std::vector< double > segmentCalls(calls.begin() + startIdx, calls.begin() + endIdx + 1);
	std::vector< double > segmentAltCalls(altCalls.begin() + startIdx, altCalls.begin() + endIdx + 1);
	
	std::vector<double> segmentDensities;
	int maxCallsInd = segmentCalls.size();
	for (int i = 0; i < maxCallsInd; i++){
	
		int positiveCalls = 0;
		int attempts = 0;
		int lb = std::max(0,i-epsilon);
		int segmentIdxLen = segmentCalls.size();
		int ub = std::min(segmentIdxLen, i+epsilon);

		for (int j = lb; j < ub; j++){
		
			int signedGap = segmentPositions[i]-segmentPositions[j];
			if (std::abs(signedGap) < epsilon){
	
				if (segmentCalls[j] > 0.5) positiveCalls++;
				if (segmentAltCalls[j] > 0.5) positiveCalls--;
				attempts++;
			}
		}
		segmentDensities.push_back((double)positiveCalls / (double) attempts);
	}
	
	
	double minDensity = vectorMean(segmentDensities);
	
	std::map<int,int> labels = DBSCAN_mod( segmentPositions, segmentCalls, segmentAltCalls, epsilon, minDensity );
	
	int trimFromLeft = 0;
	int maxPosInd = segmentPositions.size();
	for (int i = 0; i < maxPosInd; i++){
	
		if (labels.at(i) < 0) trimFromLeft++;
		else break;
	}
	
	int trimFromRight = 0;
	for (int i = maxPosInd - 1; i > 0; i--){
	
		if (labels.at(i) < 0) trimFromRight++;
		else break;
	}

	return std::make_pair(trimFromLeft, trimFromRight);
}


void callStalls(DetectedRead &r, std::string analogueOrder, KMeansResult analougeIncorporation){

	int filterSize = 2000;
	
	std::vector< double > secondAnalogueCalls;
	if (analogueOrder == "EdU,BrdU"){
	
		secondAnalogueCalls = r.brduCalls;
	}
	else{
	
		secondAnalogueCalls = r.eduCalls;
	}

	//non-linear scaling parameters for stall score
	double beta = 10.; //higher values of beta mean more conservative stall scores
	double alpha = 1./log(2./(1.+exp(-1.*beta))); //set alpha so that non-linear scaling of 1 is equal to 1

	//check right forks
	for (auto s = r.rightForks.begin(); s < r.rightForks.end(); s++){
	
		if (s -> partners > 0){
			s -> score = -1;
			continue;
		}
	
		int forkTipIdx = s -> rightmostIdx;
		int numPositions = r.positions.size();
		assert( forkTipIdx < numPositions);
		double maximumScore = -3.0;
		bool failOnAlignment = false;

		if (forkTipIdx > filterSize and forkTipIdx < numPositions-filterSize){
			
			int positiveCalls = 0;
			int attempts = 0;
			for (int j = forkTipIdx-filterSize; j < forkTipIdx; j++){
			
				if (std::abs(r.positions[forkTipIdx] - r.positions[j]) < filterSize){

					if (not r.alignmentQuality[j]) failOnAlignment = true;
			
					if (secondAnalogueCalls[j] > 0.5){
						positiveCalls++;
					}
					attempts++;
				}
			}
			if (attempts < 50) continue;
			double LHS = (double) positiveCalls / (double) attempts;
			
			//guard against low denominators
			if (LHS < 0.2) continue;
			
			positiveCalls = 0;
			attempts = 0;
			for (int j = forkTipIdx; j < forkTipIdx+filterSize; j++){
			
				if (std::abs(r.positions[forkTipIdx] - r.positions[j]) < filterSize){

					if (not r.alignmentQuality[j]) failOnAlignment = true;
			
					if (secondAnalogueCalls[j] > 0.5){
						positiveCalls++;
					}
					attempts++;
				}
			}
			if (attempts < 50) continue;
			double RHS = (double) positiveCalls / (double) attempts;
		
			double score;
			if (LHS - RHS > 0.){
				score = (LHS-RHS)/LHS;
				score = alpha*log(1+exp(beta*(score-1))) - alpha*log(1+exp(beta*(-1)));
			}
			else{
				score = -2.0;
			}
			
			r.stallScore[forkTipIdx] = score;
			if (score > maximumScore){
				maximumScore = score;
			}
		}

		if (failOnAlignment){
			s -> score = -4.0; 
			continue;
		}

		s -> score = maximumScore;
	}
	
	//check left forks
	for (auto s = r.leftForks.begin(); s < r.leftForks.end(); s++){
	
		if (s -> partners > 0){
			s -> score = -1;
			continue;
		}
	
		int forkTipIdx = s -> leftmostIdx;
		int numPositions = r.positions.size();
		assert( forkTipIdx < numPositions);
		double maximumScore = -3.0;
		bool failOnAlignment = false;

		if (forkTipIdx > filterSize and forkTipIdx < numPositions-filterSize){
			
			int positiveCalls = 0;
			int attempts = 0;
			for (int j = forkTipIdx-filterSize; j < forkTipIdx; j++){
			
				if (std::abs(r.positions[forkTipIdx] - r.positions[j]) < filterSize){

					if (not r.alignmentQuality[j]) failOnAlignment = true;
				
					if (secondAnalogueCalls[j] > 0.5){
						positiveCalls++;
					}
					attempts++;
				}
			}
			if (attempts < 50) continue;
			
			double LHS = (double) positiveCalls / (double) attempts;
			
			positiveCalls = 0;
			attempts = 0;
			for (int j = forkTipIdx; j < forkTipIdx+filterSize; j++){
			
				if (std::abs(r.positions[forkTipIdx] - r.positions[j]) < filterSize){

					if (not r.alignmentQuality[j]) failOnAlignment = true;
			
					if (secondAnalogueCalls[j] > 0.5){
						positiveCalls++;
					}
					attempts++;
				}
			}
			if (attempts < 50) continue;
			
			double RHS = (double) positiveCalls / (double) attempts;
			
			if (RHS < 0.2) continue;
		
			double score;
			if (RHS - LHS > 0.){
				score = (RHS-LHS)/RHS;
				score = alpha*log(1+exp(beta*(score-1))) - alpha*log(1+exp(beta*(-1)));
			}
			else{
				score = -2.0;
			}
			r.stallScore[forkTipIdx] = score;
			if (score > maximumScore){
				maximumScore = score;
			}
		}

		if (failOnAlignment){
			s -> score = -4.0;
			continue;
		}

		s -> score = maximumScore;
	}
}


void emptyBuffer(std::vector< DetectedRead > &buffer, forkSenseArgs args, fs_fileManager &fm, KMeansResult analogueIncorporation){

	#pragma omp parallel for schedule(dynamic) shared(args, analogueIncorporation) num_threads(args.threads)
	for ( auto b = buffer.begin(); b < buffer.end(); b++) {

		runDBSCAN(*b,analogueIncorporation, args.analogueOrder);
		callSegmentation(*b);
		
		std::string termOutput, originOutput, leftForkOutput, rightForkOutput, leftForkOutput_signatures, rightForkOutput_signatures, BrdUOutput, EdUOutput;
		bool segmentToForks = false;

		if (args.markOrigins or args.markTerms or args.markForks){

			callForks(*b, args.analogueOrder);
			callStalls(*b, args.analogueOrder, analogueIncorporation);
			
			for (auto lf = (*b).leftForks.begin(); lf < (*b).leftForks.end(); lf++){
			
				leftForkOutput += (*b).chromosome + " " 
				               + std::to_string(lf -> leftmostCoord) + " " 
				               + std::to_string(lf -> rightmostCoord) + " " 
				               + (*b).readID.substr(1) + " "
				               + std::to_string( (*b).mappingLower ) + " "
				               + std::to_string( (*b).mappingUpper ) + " "
				               + (*b).strand + " "
				               + std::to_string(lf -> score) + "\n"; 
			}
			
			for (auto rf = (*b).rightForks.begin(); rf < (*b).rightForks.end(); rf++){
			
				rightForkOutput += (*b).chromosome + " " 
				               + std::to_string(rf -> leftmostCoord) + " " 
				               + std::to_string(rf -> rightmostCoord) + " " 
				               + (*b).readID.substr(1) + " "
				               + std::to_string( (*b).mappingLower ) + " "
				               + std::to_string( (*b).mappingUpper ) + " "
				               + (*b).strand + " "
				               + std::to_string(rf -> score) + "\n"; 
			}
			
			if (args.makeSignatures){
			
				for (auto lf = (*b).leftForks.begin(); lf < (*b).leftForks.end(); lf++){
				
					leftForkOutput_signatures += (*b).chromosome + " " 
						       + std::to_string(lf -> leftmostCoord) + " " 
						       + std::to_string(lf -> rightmostCoord) + " " 
						       + (*b).readID.substr(1) + " "
						       + std::to_string( (*b).mappingLower ) + " "
						       + std::to_string( (*b).mappingUpper ) + " "
						       + (*b).strand + " ";
					for (auto s = (lf -> stress_signature).begin(); s < (lf -> stress_signature).end(); s++) leftForkOutput_signatures += std::to_string(*s) + " ";
					leftForkOutput_signatures += std::to_string(lf -> score) + "\n"; 
				}
				
				for (auto rf = (*b).rightForks.begin(); rf < (*b).rightForks.end(); rf++){
				
					rightForkOutput_signatures += (*b).chromosome + " " 
						       + std::to_string(rf -> leftmostCoord) + " " 
						       + std::to_string(rf -> rightmostCoord) + " " 
						       + (*b).readID.substr(1) + " "
						       + std::to_string( (*b).mappingLower ) + " "
						       + std::to_string( (*b).mappingUpper ) + " "
						       + (*b).strand + " ";
					for (auto s = (rf -> stress_signature).begin(); s < (rf -> stress_signature).end(); s++) rightForkOutput_signatures += std::to_string(*s) + " ";
					rightForkOutput_signatures += std::to_string(rf -> score) + "\n"; 
				}			
			}

			if (args.markOrigins){

				originOutput = callOrigins(*b,args);
			}
			if (args.markTerms){

				termOutput = callTerminations(*b,args);
			}
			segmentToForks = true;
		}
		
		if (args.markAnalogues){

			std::pair<std::string,std::string> analogueOutputPair = writeAnalogueRegions(*b, segmentToForks);
			BrdUOutput = analogueOutputPair.first;
			EdUOutput = analogueOutputPair.second;
		}
		
		//fix the analogue segment indices
		std::vector<int> eduSegment_output( (b -> positions).size(), 0);
		std::vector<int> brduSegment_output( (b -> positions).size(), 0);
		for (auto s = (*b).EdU_segment.begin(); s < (*b).EdU_segment.end(); s++){
		
			if ( (*s).partners == 0 ) continue;
		
			for (int i = (*s).leftmostIdx; i <= (*s).rightmostIdx; i++) eduSegment_output[i] = 1;
		}
		for (auto s = (*b).BrdU_segment.begin(); s < (*b).BrdU_segment.end(); s++){
		
			if ( (*s).partners == 0 ) continue;
		
			for (int i = (*s).leftmostIdx; i <= (*s).rightmostIdx; i++) brduSegment_output[i] = 1;
		}
		
		//only output segmentation on non-trivial reads that have at least one analogue segment called on them
		std::string readOutput;
		if ( (*b).EdU_segment.size() > 0 or (*b).BrdU_segment.size() > 0){
		
			//write the read header
			readOutput += (*b).readID + " " + (*b).chromosome + " " + std::to_string((*b).mappingLower) + " " + std::to_string((*b).mappingUpper) + " " + (*b).strand + "\n"; //header

			for (size_t i = 0; i < (*b).positions.size(); i++){
			
				readOutput += std::to_string((*b).positions[i]) + "\t" + std::to_string(eduSegment_output[i]) + "\t" + std::to_string(brduSegment_output[i]) + "\n";
			}
		}

		#pragma omp critical
		{
			fm.writeOutput(readOutput, termOutput, originOutput, leftForkOutput, rightForkOutput, leftForkOutput_signatures, rightForkOutput_signatures, BrdUOutput, EdUOutput);
		}
	}
	buffer.clear();
}


KMeansResult twoMeans_fs( std::vector< double > &observations ){

	double C1_old = 0.01;
	double C2_old = 0.5;
	double C1_new = C1_old;
	double C2_new = C2_old;
	double tol = 0.0001;
	int maxIter = 100;
	int iter = 0;

	std::vector<double> C1_points_old;
	std::vector<double> C2_points_old;

	//make an initial assignment
	for ( size_t i = 0; i < observations.size(); i++ ){

		if ( std::abs(observations[i] - C1_old) < std::abs(observations[i] - C2_old) ) C1_points_old.push_back(observations[i]);
		else C2_points_old.push_back(observations[i]);
	}

	//iterate until tolerance is met
	do{
		C1_old = C1_new;
		C2_old = C2_new;

		std::vector<double> C1_points_new;
		std::vector<double> C2_points_new;

		for ( size_t i = 0; i < C1_points_old.size(); i++ ){

			if ( std::abs(C1_points_old[i] - C1_old) < std::abs(C1_points_old[i] - C2_old) ) C1_points_new.push_back(C1_points_old[i]);
			else C2_points_new.push_back(C1_points_old[i]);
		}	

		for ( size_t i = 0; i < C2_points_old.size(); i++ ){

			if ( std::abs(C2_points_old[i] - C1_old) < std::abs(C2_points_old[i] - C2_old) ) C1_points_new.push_back(C2_points_old[i]);
			else C2_points_new.push_back(C2_points_old[i]);
		}	

		C1_new = vectorMean(C1_points_new);
		C2_new = vectorMean(C2_points_new);

		C1_points_old = C1_points_new;
		C2_points_old = C2_points_new;

		iter++;
	}while (iter < maxIter and (std::abs(C1_old - C1_new)>tol or std::abs(C2_old - C2_new)>tol));
	
	//calculate lower bounds from K-means segmentation
	auto C1_lowerBound = std::min_element( C1_points_old.begin(), C1_points_old.end() );
	auto C2_lowerBound = std::min_element( C2_points_old.begin(), C2_points_old.end() );

	//more conservative option - lower bound is mean - 1 stdv
	double C1_stdv = vectorStdv( C1_points_old, C1_new );
	double C2_stdv = vectorStdv( C2_points_old, C2_new );
	
	KMeansResult out = {C1_new, *C1_lowerBound, C1_stdv, C2_new, *C2_lowerBound, C2_stdv};

	return out;
}


int parseDetectLine(std::string line,
			int &BrdUcalls,
			int &EdUcalls,
			int &attempts){

	std::string column;
	std::stringstream ssLine(line);
	int position = -1, cIndex = 0;
	AnalogueScore B, E;
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
		cIndex++;
	}
	assert(position != -1);


	if ( B.get() > 0.5 ){
		BrdUcalls++;
		attempts++;
	}
	else if ( E.get() > 0.5 ){
		EdUcalls++;
		attempts++;
	}
	else{
		attempts++;
	}
	return position;
}

KMeansResult estimateAnalogueIncorporation(std::string detectFilename, int readCount){

	int readCap = readCount;

	std::ifstream inFile( detectFilename );
	if ( not inFile.is_open() ) throw IOerror( detectFilename );

	std::cout << "Estimating analogue incorporation..." << std::endl;

	std::vector< double > BrdU_callFractions, EdU_callFractions;
	
	int resolution = 1000; //look in 1 kb segments
	
	int startingPos = -1;
	int progress = 0;
	std::string line;
	
	progressBar pb(readCap,false);
	
	int BrdUcalls = 0;
	int EdUcalls = 0;
	int attempts = 0;
	int gap = 0;
	
	while( std::getline( inFile, line ) ){

		if (line.substr(0,1) == "#") continue; //ignore header
		if ( line.substr(0,1) == ">" ){

			progress++;
			pb.displayProgress( progress, 0, 0 );
			
			if (progress >= readCap) break;
			
			BrdUcalls = 0, EdUcalls = 0, attempts = 0, gap = 0, startingPos = -1;
			continue;
		}

		int position = parseDetectLine(line, BrdUcalls, EdUcalls, attempts);

		if (position == -1) continue;

		if ( startingPos == -1 ) startingPos = position;
		gap = position - startingPos;

		if ( gap > resolution and attempts >= resolution / 10 ){

			double BrdU_frac = (double) BrdUcalls / (double) attempts;
			BrdU_callFractions.push_back( BrdU_frac );

			double EdU_frac = (double) EdUcalls / (double) attempts;
			EdU_callFractions.push_back( EdU_frac );

			BrdUcalls = 0, EdUcalls = 0, attempts = 0, gap = 0, startingPos = -1;
		}
	}
	std::cout << std::endl << "Done." << std::endl;

	KMeansResult BrdU_KMeans = twoMeans_fs( BrdU_callFractions );
	
	double BrdU_p;
	double BrdU_stdv;
	double BrdU_lowerBound;
	
	if (BrdU_KMeans.centroid_1 > BrdU_KMeans.centroid_2){
		BrdU_p = BrdU_KMeans.centroid_1;
		BrdU_stdv = BrdU_KMeans.centroid_1_stdv;
		BrdU_lowerBound = BrdU_KMeans.centroid_1_lowerBound;
	}
	else{
		BrdU_p = BrdU_KMeans.centroid_2;
		BrdU_stdv = BrdU_KMeans.centroid_2_stdv;
		BrdU_lowerBound = BrdU_KMeans.centroid_2_lowerBound;	
	}
	
	
	KMeansResult EdU_KMeans = twoMeans_fs( EdU_callFractions );
	
	double EdU_p;
	double EdU_stdv;
	double EdU_lowerBound;
	
	if (EdU_KMeans.centroid_1 > EdU_KMeans.centroid_2){
		EdU_p = EdU_KMeans.centroid_1;
		EdU_stdv = EdU_KMeans.centroid_1_stdv;
		EdU_lowerBound = EdU_KMeans.centroid_1_lowerBound;
	}
	else{
		EdU_p = EdU_KMeans.centroid_2;
		EdU_stdv = EdU_KMeans.centroid_2_stdv;
		EdU_lowerBound = EdU_KMeans.centroid_2_lowerBound;	
	}

	std::cerr << "Estimated fraction of BrdU substitution in BrdU-positive regions: " << BrdU_p << std::endl;
	std::cerr << "Estimated BrdU substitution lower bound in BrdU-positive regions: " << BrdU_lowerBound << std::endl;
	std::cerr << "Estimated fraction of EdU substitution in EdU-positive regions: " << EdU_p << std::endl;
	std::cerr << "Estimated EdU substitution lower bound in EdU-positive regions: " << EdU_lowerBound << std::endl;

	inFile.close();
	
	KMeansResult out = {BrdU_p, BrdU_lowerBound, BrdU_stdv, EdU_p, EdU_lowerBound, EdU_stdv};
	
	return out;
}


int sense_main( int argc, char** argv ){

	forkSenseArgs args = parseSenseArguments( argc, argv );

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
 	inFile.open( args.detectFilename );
	if ( not inFile.is_open() ) throw IOerror( args.detectFilename );
 	std::ofstream outFile( args.outputFilename );
 	
 	fs_fileManager fm(args, analogueIncorporation);

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
			if (readBuffer.size() >= maxBufferSize) emptyBuffer(readBuffer, args, fm, analogueIncorporation);

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
				else if ( cIndex == 4 ){
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
	emptyBuffer(readBuffer, args, fm, analogueIncorporation);

	inFile.close();
	fm.closeAll();

	std::cout << std::endl << "Done." << std::endl;

	return 0;
}

