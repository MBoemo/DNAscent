//----------------------------------------------------------
// Copyright 2025 University of Cambridge
// This software is licensed under GPL-3.0.  You should have
// received a copy of the license with this software.  If
// not, please Email the author.
//----------------------------------------------------------

#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#include <unordered_set>
#include <cstdlib>  
#include <ctime>    
#include <cmath>
#include "seeBreaks.h"
#include "common.h"
#include "error_handling.h"
#include "htsInterface.h"
#include "event_handling.h"
#include <numeric>
#include <random>

static const char *help=
"seeBreaks: DNAscent executable that detects an elevated frequency of DNA breaks at replication forks.\n"
"To run DNAscent seeBreaks, do:\n"
"   DNAscent seeBreaks -l /path/to/leftForks_DNAscent_forksense.bed -r /rightForks_DNAscent_forksense.bed -d /path/to/detectOutput.bam -o /path/to/output.seeBreaks \n"
"Required arguments are:\n"
"  -l,--left                 path to leftForks file from forkSense detect with `bed` extension,\n"
"  -r,--right                path to rightFork file from forkSense detect with `bed` extension,\n"
"  -a,--analogue             path to second pulsed analogue file from forkSense with `bed` extension,\n"
"  -d,--detect               path to output from detect with `detect` or `bam` extension,\n"
"  -o,--output               path to output file or directory for seeBreaks,\n"
"DNAscent is under active development by the Boemo Group, Department of Pathology, University of Cambridge (https://www.boemogroup.org/).\n"
"Please submit bug reports to GitHub Issues (https://github.com/MBoemo/DNAscent/issues).";


struct Arguments {
	std::string lForkInput;
	std::string rForkInput;
	std::string analogueInput;
	std::string DetectInput;
	std::string output;
    bool specifiedLeft = false;
    bool specifiedRight = false;
    bool specifiedAnalogue = false;
    bool humanReadable = false;
	bool specifiedOutput = false;
    bool specifiedDetect = false;
};


Arguments parseBreaksArguments( int argc, char** argv ) {

    if( argc < 2 ){
 		std::cout << "Exiting with error.  Insufficient arguments passed to DNAscent seeBreaks." << std::endl << help << std::endl;
		exit(EXIT_FAILURE);
	}
 	if ( std::string( argv[ 1 ] ) == "-h" or std::string( argv[ 1 ] ) == "--help" ){
 		std::cout << help << std::endl;
		exit(EXIT_SUCCESS);
	}

	Arguments args;

	for ( int i = 1; i < argc; ){

 		std::string flag( argv[ i ] );

 		if ( flag == "-l" or flag == "--left" ){

 			if (i == argc-1) throw TrailingFlag(flag);

            std::string strArg( argv[ i + 1 ] );
            const char *ext = get_ext(strArg.c_str());

            if (strcmp(ext,"bed") != 0) throw InvalidExtension(ext);

			args.lForkInput = strArg;
			i+=2;
			args.specifiedLeft = true;
		}
        else if ( flag == "-r" or flag == "--right" ){

 			if (i == argc-1) throw TrailingFlag(flag);
            
            std::string strArg( argv[ i + 1 ] );
            const char *ext = get_ext(strArg.c_str());
            
            if (strcmp(ext,"bed") != 0) throw InvalidExtension(ext);
 					
			args.rForkInput = strArg;
			i+=2;
			args.specifiedRight = true;
        }
        else if ( flag == "-a" or flag == "--analogue" ){

 			if (i == argc-1) throw TrailingFlag(flag);
            
            std::string strArg( argv[ i + 1 ] );
            const char *ext = get_ext(strArg.c_str());
            
            if (strcmp(ext,"bed") != 0) throw InvalidExtension(ext);
 					
			args.analogueInput = strArg;
			i+=2;
			args.specifiedAnalogue = true;
        }
        else if ( flag == "-h" or flag == "--help" ){
            std::cout << help << std::endl;
            exit(EXIT_SUCCESS);
		}
        else if ( flag == "-d" or flag == "--detect" ){
 		
 			if (i == argc-1) throw TrailingFlag(flag);		
 		
 			std::string strArg( argv[ i + 1 ] );
 			
 			const char *ext = get_ext(strArg.c_str());
				
			if (strcmp(ext,"bam") == 0){
				args.humanReadable = false;
			}
			else if (strcmp(ext,"detect") == 0){
				args.humanReadable = true;
            }
			else{
				throw InvalidExtension(ext);
			}
            args.specifiedDetect = true;
			
			args.DetectInput = strArg;
			i+=2;
		}
		else if ( flag == "-o" or flag == "--output" ){
		
			if (i == argc-1) throw TrailingFlag(flag);		
		
 			std::string strArg( argv[ i + 1 ] );
            args.output = strArg;
			i+=2;
			args.specifiedOutput = true;
        }
        else throw InvalidOption( flag );
	}
    if ( !args.specifiedOutput or !args.specifiedAnalogue or !args.specifiedDetect or (!args.specifiedLeft and !args.specifiedRight) ) {
        std::cout << "Exiting with error.  Insufficient arguments passed to DNAscent seeBreaks." << std::endl;
        exit(EXIT_FAILURE);
    }

    return args;
}


void detectUnpack(Arguments &args, std::vector<int> &v5Prime, std::vector<int> &v3Prime, int minReadLength) {

    std::ifstream detectFile(args.DetectInput);
    std::string line;
    int prog = 0;

	while( std::getline( detectFile, line ) ){

        if (line.substr(0,1) == "#" or line.length() == 0) continue; //ignore header and blank lines
        if ( line.substr(0,1) == ">" ){
		
            //parse the header line
            std::istringstream ssLine(line);
            std::string column;
            int cIndex = 0;
            int refStart = -1;
            int refEnd = -1;
            while ( ssLine >> column ) {

                if ( cIndex == 2 ) {

                    refStart = std::stoi(column);
                }
                else if ( cIndex == 3 ) {

                    refEnd = std::stoi(column);
                }
                cIndex++;
            }
            assert(refStart != -1 and refEnd != -1);

            if (refEnd - refStart < minReadLength) continue; //ignore short reads
            
            v5Prime.push_back(refStart);
            v3Prime.push_back(refEnd);

            prog ++;
            if (prog % 1000 == 0) std::cout << "\rProcessed " << prog << " reads..." << std::flush;	
        }
    }
    detectFile.close();
}


void bamUnpack(Arguments &args, std::vector<int> &v5Prime, std::vector<int> &v3Prime, int minReadLength) {

    htsFile *bam_fh = sam_open((args.DetectInput).c_str(), "r");
    if (bam_fh == NULL) throw IOerror(args.DetectInput);

    bam_hdr_t *bam_hdr = sam_hdr_read(bam_fh);
    bam1_t *itr_record = bam_init1();

    int prog = 0;

    while(sam_read1(bam_fh, bam_hdr, itr_record) >= 0){
          
        int refStart,refEnd;
        getRefEnd(itr_record,refStart,refEnd);

        if (refEnd - refStart < minReadLength) continue; //ignore short reads

        v5Prime.push_back(refStart);
        v3Prime.push_back(refEnd);

        prog ++;
        if (prog % 1000 == 0) std::cout << "\rProcessed " << prog << " reads..." << std::flush;	
    }
	bam_destroy1(itr_record);				
	bam_hdr_destroy(bam_hdr);
	hts_close(bam_fh);
}


void scanReadIDs(std::string fileInput, std::vector<std::string> &ReadIDs, std::vector<std::string> &DuplicateIDs ) {
    
    std::ifstream file(fileInput);

    if (!file.is_open()) 
    {
        std::cerr << "Error: Could not open file " << fileInput << "\n";
        exit(EXIT_FAILURE);
    }
    
    std::string line; 

    while (std::getline(file, line)) {

        if (!line.empty() && line[0] != '#') {

            std::istringstream iss(line);
            std::vector<std::string> columns;
            std::string entry;

            while (iss >> entry) {

                columns.push_back(entry);
            }
                                
            std::string readID = columns[3];
            if (std::find(ReadIDs.begin(), ReadIDs.end(), readID) != ReadIDs.end()) {

                DuplicateIDs.push_back(readID);
            }
            else{

                ReadIDs.push_back(readID);
            }
        }
    }
    file.close();
}


void analogueUnpack(std::string fileInput, std::map<std::string, std::pair<int,int>> &readID2analogue, std::vector<std::string> &DuplicateIDs) {
    
    std::ifstream file(fileInput);

    if (!file.is_open()) 
    {
        std::cerr << "Error: Could not open file " << fileInput << "\n";
        exit(EXIT_FAILURE);
    }
    
    std::string line; 

    while (std::getline(file, line)) {

        if (!line.empty() && line[0] != '#') {

            std::istringstream iss(line);
            std::vector<std::string> columns;
            std::string entry;
            while (iss >> entry) {
                columns.push_back(entry);
            }
                                
            std::string readID = columns[3];

            if ( std::find(DuplicateIDs.begin(), DuplicateIDs.end(), readID) == DuplicateIDs.end() ) {

                int pulse5Prime = std::stoi(columns[1]);
                int pulse3Prime = std::stoi(columns[2]);

                readID2analogue[readID] = std::make_pair(pulse5Prime, pulse3Prime);
            }
        }
    }
    file.close();
}


void getForkSpeed(std::string fileInput, bool isRight, std::vector<int> &ForkLength, std::map<std::string, std::pair<int,int>> &readID2analogue, std::vector<std::string> &duplicateIDs, int &fsBoundary, bool &isR9) {
    
    std::ifstream file(fileInput);

    if (!file.is_open()) 
    {
        std::cerr << "Error: Could not open file " << fileInput << "\n";
        exit(EXIT_FAILURE);
    }
    
    std::string line; 

    while (std::getline(file, line)) {

        if (!line.empty() && line[0] != '#') {

            std::istringstream iss(line);
            std::vector<std::string> columns;
            std::string entry;
            while (iss >> entry) {
                columns.push_back(entry);
            }

            if (columns.size() == 8) isR9 = true;
            else if (columns.size() == 9) isR9 = false;
            else{
                std::cerr << "Error: Incorrect number of columns in forkSense bed file " << fileInput << "\n";
                exit(EXIT_FAILURE);
            }
                                
            std::string readID = columns[3];

            if ( std::find(duplicateIDs.begin(), duplicateIDs.end(), readID) == duplicateIDs.end() ) {

                int pulse5Prime = std::stoi(columns[1]);
                int pulse3Prime = std::stoi(columns[2]);
                int read5Prime = std::stoi(columns[4]);
                int read3Prime = std::stoi(columns[5]);

                int secondAn3Prime;
                int secondAn5Prime;

                if (isRight and pulse3Prime == readID2analogue[readID].second){
                   secondAn3Prime = readID2analogue[readID].second;
                   secondAn5Prime = readID2analogue[readID].first;
                }
                else if (not isRight and pulse5Prime == readID2analogue[readID].first){
                   secondAn3Prime = readID2analogue[readID].second;
                   secondAn5Prime = readID2analogue[readID].first;
                }
                else{
                    continue;
                }

                int gap3Prime = read3Prime - secondAn3Prime;
                int gap5Prime = secondAn5Prime - read5Prime;
                assert(gap3Prime >= 0 and gap5Prime >= 0);

                // For reliable fork speeds, only consider forks that are at least fsBoundary bp away from the read ends
                if ( (gap3Prime > fsBoundary) and (gap5Prime > fsBoundary ) ) {

                    ForkLength.push_back(secondAn3Prime - secondAn5Prime);
                }
            }
        }
    }
    file.close();
}


int forkUnpack(std::string fileInput, bool isRight, std::map<std::string, std::pair<int,int>> &readID2analogue, std::vector<bool> &runOff, std::vector<std::string> &duplicateIDs, int &fsBoundary, int &readEndTolerance, int &minReadLength) {
    
    std::ifstream file(fileInput);
    int nForks = 0;

    if (!file.is_open()) 
    {
        std::cerr << "Error: Could not open file " << fileInput << "\n";
        exit(EXIT_FAILURE);
    }
    
    std::string line; 

    while (std::getline(file, line)) {

        if (!line.empty() && line[0] != '#') {

            std::istringstream iss(line);
            std::vector<std::string> columns;
            std::string entry;
            while (iss >> entry) {
                columns.push_back(entry);
            }
                                
            std::string readID = columns[3];

            if ( std::find(duplicateIDs.begin(), duplicateIDs.end(), readID) == duplicateIDs.end() ) {

                int pulse5Prime = std::stoi(columns[1]);
                int pulse3Prime = std::stoi(columns[2]);
                int read5Prime = std::stoi(columns[4]);
                int read3Prime = std::stoi(columns[5]);

                if (read3Prime - read5Prime < minReadLength) continue; //ignore short reads

                int secondAn3Prime;
                int secondAn5Prime;

                if (isRight and pulse3Prime == readID2analogue[readID].second){
                   secondAn3Prime = readID2analogue[readID].second;
                   secondAn5Prime = readID2analogue[readID].first;
                }
                else if (not isRight and pulse5Prime == readID2analogue[readID].first){
                   secondAn3Prime = readID2analogue[readID].second;
                   secondAn5Prime = readID2analogue[readID].first;
                }
                else{
                    continue;
                }

                int gap3Prime = read3Prime - secondAn3Prime;
                int gap5Prime = secondAn5Prime - read5Prime;
                assert(gap3Prime >= 0 and gap5Prime >= 0);

                if ( isRight and (gap5Prime > fsBoundary)){
                    if (gap3Prime < readEndTolerance) runOff.push_back(true);
                    else runOff.push_back(false);
                    nForks++;  
                }
                else if ( not isRight and (gap3Prime > fsBoundary) ){
                    if (gap5Prime < readEndTolerance) runOff.push_back(true);
                    else runOff.push_back(false);
                    nForks++;               
                }
            }
        }
    }
    file.close();
    return nForks;
}


void simulation (std::vector<int> &v5Prime, 
                std::vector<int> &v3Prime, 
                std::vector<int> &forkLength,
                unsigned int nForks,
                std::vector<double> &totalRunOffs,
                int fsBoundary,
                int readEndTolerance) { 

    int bsIterations = 5000;
    std::mt19937 gen(221005);

    // Fork simulation
    for (int i = 0; i < bsIterations; i++) {  

        int runOff = 0;           
            
        for (size_t j = 0; j < nForks; ++j) {

            // Select boundaries of a random read
            assert(v5Prime.size() == v3Prime.size());
            std::uniform_int_distribution<> readDist(0, v5Prime.size()-1);
            int readIndex = readDist(gen);
            int read5Prime = v5Prime[readIndex];
            int read3Prime = v3Prime[readIndex];

            // Randomly select a fork track length
            std::uniform_int_distribution<> trackDist(0 , forkLength.size()-1);
            int trackIndex = trackDist(gen);
            int randomLength = forkLength[trackIndex];

            // Randomly select a starting point on the read
            int upperCutoff = read3Prime - fsBoundary;
            int lowerCutoff = read5Prime + fsBoundary;
            std::uniform_int_distribution<> startDist(lowerCutoff, upperCutoff);
            int randomStart = startDist(gen);

            // Check if there is a run off
            if (read3Prime - readEndTolerance - randomStart < randomLength) runOff++;  
        }

        // Convert to proportion and store proportion
        double probRunOff = static_cast<double>(runOff) / nForks;
        totalRunOffs.push_back(probRunOff);
    }
}

void observation(std::vector<bool> &runOff, std::vector<double> &totalObsRunOff) {

    int bsIterations = 5000;
    std::mt19937 gen(221005);

    for (int i = 0; i < bsIterations; ++i) {

        int obsRunOffs = 0;
        int noObsRunOffs = 0;

        for (size_t j = 0; j < runOff.size(); ++j) {

            std::uniform_int_distribution<> distrib(0 , runOff.size()-1);
            int randomIndex = distrib(gen);
            bool randomRunOff = runOff[randomIndex];

            if (randomRunOff) {
                obsRunOffs += 1;
            }
            else {
                noObsRunOffs += 1;
            }
        }
        double probObsRunOff = static_cast<double>(obsRunOffs) / (obsRunOffs + noObsRunOffs);
        totalObsRunOff.push_back(probObsRunOff);
    }
}


int seeBreaks_main(int argc, char** argv) {

    Arguments args = parseBreaksArguments(argc, argv);

    int forkSenseBoundary = 2000;
    bool isR9 = false;

    std::vector<std::string> ReadIDs;
    std::vector<std::string> DuplicateIDs;
    
    if (args.specifiedLeft) {

        scanReadIDs(args.lForkInput, ReadIDs, DuplicateIDs);
    }
    if (args.specifiedRight) {

        scanReadIDs(args.rForkInput, ReadIDs, DuplicateIDs);
    }
    
    std::map<std::string, std::pair<int,int>> readID2analogue;
    analogueUnpack(args.analogueInput, readID2analogue, DuplicateIDs);

    std::vector<int> forkTrackLengths;
    if (args.specifiedLeft) {

        getForkSpeed(args.lForkInput, false, forkTrackLengths, readID2analogue, DuplicateIDs, forkSenseBoundary, isR9);
    }
    if (args.specifiedRight) {

        getForkSpeed(args.rForkInput, true, forkTrackLengths, readID2analogue, DuplicateIDs, forkSenseBoundary, isR9);
    }

    double meanForkSpeed = vectorMean(forkTrackLengths);
    double stdvForkSpeed = vectorStdv(forkTrackLengths, meanForkSpeed);
    int minReadLength = (int) meanForkSpeed + 3 * stdvForkSpeed;

    // Parse detect output
    std::vector<int> v5Prime;
    std::vector<int> v3Prime;

    if (args.humanReadable) {

        detectUnpack(args, v5Prime, v3Prime, minReadLength);
    }
    else{// is modbam
        
        bamUnpack(args, v5Prime, v3Prime, minReadLength);       
    }

    // Parse forkSense bed files
    std::vector<double> totalSimRunOffs;
    std::vector<double> totalObsRunOffs;
    int nRightForks = 0;
    int nLeftForks = 0;

    int endTol = 250; //R10
    if (isR9) endTol = 500; //adjust up for higher epsilon in R9 v3.1.2 DBSCAN

    for (int readEndTolerance = endTol; readEndTolerance <= endTol + 250; readEndTolerance += 50) {

        std::vector<bool> runOffs;

        if (args.specifiedLeft) {

            nLeftForks = forkUnpack(args.lForkInput, false, readID2analogue, runOffs, DuplicateIDs, forkSenseBoundary, readEndTolerance, minReadLength);
        }
        if (args.specifiedRight) {

            nRightForks = forkUnpack(args.rForkInput, true, readID2analogue, runOffs, DuplicateIDs, forkSenseBoundary, readEndTolerance, minReadLength);
        }

        simulation(v5Prime, v3Prime, forkTrackLengths, runOffs.size(), totalSimRunOffs, forkSenseBoundary, readEndTolerance);
        observation(runOffs, totalObsRunOffs);  
    }
    int nForks = nRightForks + nLeftForks;

    // Calculate simulation mean and standard deviation
    double simMean = vectorMean(totalSimRunOffs);
    double simStdDev = vectorStdv(totalSimRunOffs, simMean);

    // Calculate observed mean and standard deviation  
    double obsMean = vectorMean(totalObsRunOffs);
    double obsStdDev = vectorStdv(totalObsRunOffs, obsMean);

    // Difference between simulated and observed
    std::mt19937 gen(221005);
    std::vector<double> difference;

    for (size_t i = 0; i < totalSimRunOffs.size(); ++i) {
        std::normal_distribution<double> obsDistribution(obsMean, obsStdDev);
        std::normal_distribution<double> simDistribution(simMean, simStdDev);   
        difference.push_back(obsDistribution(gen) - simDistribution(gen));
    }
    double difMean = vectorMean(difference);
    double difStdDev = vectorStdv(difference, difMean);
    double leftTail = difMean - 1.96 * difStdDev;
    double rightTail = difMean + 1.96 * difStdDev;

    // Write output to stdout
    std::cout << "\nNumber of forks: " << nForks << "\n";
    std::cout << "Expected number of analogue tracks at read ends\n";
    std::cout << "   Estimate: " << simMean << "\n";
    std::cout << "   StandardError: " << simStdDev << "\n";
    std::cout << "Observed number of analogue tracks at read ends\n";
    std::cout << "   Estimate: " << obsMean << "\n";
    std::cout << "   StandardError: " << obsStdDev << "\n";
    std::cout << "Difference between observed and expected\n";
    std::cout << "   Estimate: " << difMean << "\n"; 
    std::cout << "   StandardError: " << difStdDev << "\n";
    std::cout << "   95% Confidence Interval: [" << leftTail << ", " << rightTail << "]\n";

    // Write output to file
    auto t = std::time(nullptr);
	auto tm = *std::localtime(&t);
	std::ostringstream oss;
	oss << std::put_time(&tm, "%d/%m/%Y %H:%M:%S");
	auto startTime = oss.str();
    std::ofstream outFile(args.output);
    outFile << "#DetectFile " << args.DetectInput << "\n";
    outFile << "#ForkFiles " << args.lForkInput << " " << args.rForkInput << "\n";
	outFile << "#SystemStartTime " + startTime + "\n";
    outFile << "#Software " << std::string(getExePath()) << "\n";
    outFile << "#Version " << std::string(VERSION) << "\n";
    outFile << "#Commit " << std::string(getGitCommit()) << "\n";
    outFile << "#nForks " << nForks << "\n";
    outFile << "#ExpectedReadEndFraction " << simMean << "\n";
    outFile << "#ExpectedReadEndFraction_StdErr " << simStdDev << "\n";
    outFile << "#ObservedReadEndFraction " << obsMean << "\n";
    outFile << "#ObservedReadEndFraction_StdErr " << obsStdDev << "\n";
    outFile << "#Difference " << difMean << "\n";
    outFile << "#Difference_StdErr " << difStdDev << "\n";
    outFile << "#95ConfidenceInterval " << leftTail << " " << rightTail << "\n";
    outFile << ">ExpectedReadEndFractions:\n";
    for (const auto& val : totalSimRunOffs) {

        outFile << val << "\n";
    }
    outFile << ">ObservedReadEndFractions:\n";
    for (const auto& val : totalObsRunOffs) {

        outFile << val << "\n";
    }
    outFile.close();

    return 0;
}
