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
"  -d,--detect               path to output from detect with `detect` or `bam` extension,\n"
"  -o,--output               path to output file or directory for seeBreaks,\n"
"DNAscent is under active development by the Boemo Group, Department of Pathology, University of Cambridge (https://www.boemogroup.org/).\n"
"Please submit bug reports to GitHub Issues (https://github.com/MBoemo/DNAscent/issues).";


struct Arguments {
	std::string lForkInput;
	std::string rForkInput;
	std::string DetectInput;
	std::string output;
    bool specifiedLeft = false;
    bool specifiedRight = false;
    bool specifiedDetect = false;
    bool specifiedBam = false;
	bool specifiedOutput = false; 
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

    int lPresent = 0;
    int rPresent = 0;

	Arguments args;

	for ( int i = 1; i < argc; ){

 		std::string flag( argv[ i ] );

 		if ( flag == "-l" or flag == "--left" ){

            lPresent++;

 			std::string strArg( argv[ i + 1 ] );
 			const char *ext = get_ext(strArg.c_str());

            if (lPresent > 1) throw InvalidExtension(ext);
 			if (i == argc-1) throw TrailingFlag(flag);		
 		
			args.lForkInput = strArg;
			i+=2;
			args.specifiedLeft = true;
		}
        else if ( flag == "-r" or flag == "--right" ){

            rPresent++;

 			std::string strArg( argv[ i + 1 ] );
            const char *ext = get_ext(strArg.c_str());

            if (rPresent > 1) throw InvalidExtension(ext);
 			if (i == argc-1) throw TrailingFlag(flag);		
 					
			args.rForkInput = strArg;
			i+=2;
			args.specifiedRight = true;
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
				args.specifiedBam = true;
			}
			else if (strcmp(ext,"detect") == 0){
				args.specifiedDetect = true;
            }
            else if (strcmp(ext,"forkSense") == 0){
                args.specifiedDetect = true;
			}
			else{
				throw InvalidExtension(ext);
			}
			
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
    if ( !args.specifiedOutput or (!args.specifiedLeft and !args.specifiedRight) or (!args.specifiedDetect and !args.specifiedBam) ) {
        std::cout << "Exiting with error.  Insufficient arguments passed to DNAscent seeBreaks." << std::endl;
        exit(EXIT_FAILURE);
    }

    return args;
}


void detectUnpack(Arguments &args, std::vector<int> &v5Prime, std::vector<int> &v3Prime, int &dnCounter) {

    int minReadLength = 3000;
    std::ifstream detectFile(args.DetectInput);
    std::string line;

	while( std::getline( detectFile, line ) ){

        if (line.substr(0,1) == "#" or line.length() == 0) continue; //ignore header and blank lines
        if ( line.substr(0,1) == ">" ){

            dnCounter++;			
            
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
        }
    }
    detectFile.close();
}


void bamUnpack (Arguments &args, std::vector<int> &v5Prime, std::vector<int> &v3Prime, int &dnCounter) {

    int minReadLength = 3000;

    htsFile *bam_fh = sam_open((args.DetectInput).c_str(), "r");

    if (bam_fh == NULL) throw IOerror(args.DetectInput);

    bam_hdr_t *bam_hdr = sam_hdr_read(bam_fh);
    bam1_t *itr_record = bam_init1();

    while(sam_read1(bam_fh, bam_hdr, itr_record) >= 0){
                        
        int refStart,refEnd;
        getRefEnd(itr_record,refStart,refEnd);

        if (refEnd - refStart < minReadLength) continue; //ignore short reads

        v5Prime.push_back(refStart);
        v3Prime.push_back(refEnd);
    }
	bam_destroy1(itr_record);				
	bam_hdr_destroy(bam_hdr);
	hts_close(bam_fh);
}


void forkUnpack(std::string input, Arguments &args, std::vector<int> &ForkLength, std::vector<double> &StallScore, int &nCounter) {

    std::string fileInput = (input == "left") ? args.lForkInput : args.rForkInput;
    
    std::ifstream file(fileInput);

    if (!file.is_open()) 
    {
        std::cerr << "Error: Could not open file " << fileInput << "\n";
        exit(EXIT_FAILURE);
    }
    
    std::unordered_set<std::string> ReadIDs; 
    std::string line; 

    while (std::getline(file, line)) {

        if (!line.empty() && line[0] != '#') {

            std::istringstream iss(line);
            std::vector<std::string> parts;
            std::string word;

            while (iss >> word) {
                parts.push_back(word);
            }
                                
            double stallScore = 0.0;
            std::string readID = parts[3];

            if (ReadIDs.find(readID) == ReadIDs.end()) {

                int pulse5Prime = std::stoi(parts[1]);
                int pulse3Prime = std::stoi(parts[2]);

                int pulseLength = pulse3Prime - pulse5Prime;
                                    
                if (parts.size() == 8) {  // R9

                    stallScore = std::stod(parts[7]);
                } 

                else if (parts.size() == 9)  { //R10

                    stallScore = std::stod(parts[8]);
                }
                else {
                    std::cerr << "Error: Unexpected number of columns in the input file." << std::endl;
                    continue;
                }

                if (stallScore != -3) {

                    ForkLength.push_back(pulseLength);
                }

                StallScore.push_back(stallScore);
                ReadIDs.insert(readID);
                nCounter++; 
            }
        }
    }
    file.close();
}


void simulation (std::vector<int> &v5Prime, 
                std::vector<int> &v3Prime, 
                std::vector<int> &forkLength,
                std::vector<double> &stallScore,
                bool specifiedDetect, 
                bool specifiedBam,
                bool runLeft, 
                bool runRight,
                std::vector<double> &totalRunOffs) { 

    // Random number generator

    std::mt19937 gen(221005);

    // Fork simulation
    for (int i = 0; i < 5000; i++) {  
            
        int runOff = 0;
            
        for (size_t j = 0; j < stallScore.size(); ++j) {
                
            int read5Prime = 0;
            int read3Prime = 0;

            std::uniform_int_distribution<> distrib(0, v5Prime.size()-1);
            int randomIndex = distrib(gen);
            read5Prime = v5Prime[randomIndex];
            read3Prime = v3Prime[randomIndex];

            if (runLeft == true) {

                // Randomly select a starting point on the read
                std::uniform_int_distribution<> distrib(read5Prime + 2000 , read3Prime-1);
                int randomStart = distrib(gen);

                // Randomly select a pulse length
                std::uniform_int_distribution<> random(0 , forkLength.size()-1);
                int randomIndex2 = random(gen);
                int randomLength = forkLength[randomIndex2];

                // Check if there is a run off
                if (randomStart - read5Prime < randomLength) runOff++;
            }
            else if (runRight == true) {

                // Randomly select a starting point on the read
                std::uniform_int_distribution<> distrib(read5Prime , read3Prime -2000);
                int randomStart = distrib(gen);

                // Randomly select a pulse length 
                std::uniform_int_distribution<> random(1 , forkLength.size()-1);
                int randomIndex2 = random(gen);
                int randomLength = forkLength[randomIndex2];

                // Check if there is a run off
                if (read3Prime - randomStart < randomLength) runOff++;  
            }
        }

        // Convert to proportion and store proportion
        double probRunOff = static_cast<double>(runOff) / stallScore.size();
        totalRunOffs.push_back(probRunOff);
    }
}

void observation(std::vector<double> stallScore, std::vector<double> &totalObsRunOff) {

    std::random_device rd;
    std::mt19937 gen(221005);

    for (int i = 0; i < 5000; ++i) {

        int obsRunOffs = 0;
        int noObsRunOffs = 0;

        for (size_t j = 0; j < (trunc(stallScore.size() / 4)); ++j) {

            std::uniform_int_distribution<> distrib(0 , stallScore.size()-1);
            int randomIndex = distrib(gen);
            double randomScore = stallScore[randomIndex];

            if (randomScore == -3.0) {
                obsRunOffs += 1;
            }
            else if (randomScore >= 0.0) {
                noObsRunOffs += 1;
            }
        }
        double probObsRunOff = static_cast<double>(obsRunOffs) / (obsRunOffs + noObsRunOffs);
        totalObsRunOff.push_back(probObsRunOff);
    }
}


int seeBreaks_main(int argc, char** argv) {

    Arguments args = parseBreaksArguments(argc, argv);

    // Read left forks file
    int lnCounter = 0;
    std::vector<int> lForkLength;
    std::vector<double> lStallScore;

    if (args.specifiedLeft == true) {

        std::string input = "left";
        forkUnpack("left", args, lForkLength, lStallScore, lnCounter);
    }

    // Read right forks file
    std::vector<int> rForkLength;
    std::vector<double> rStallScore;
    int rnCounter = 0;

    if (args.specifiedRight == true) {

        std::string input = "right";
        forkUnpack("right", args, rForkLength, rStallScore, rnCounter);
    }

    // Parse detect output
    std::vector<int> v5Prime;
    std::vector<int> v3Prime;
    int dnCounter = 0;

    if (args.specifiedDetect == true) {

        detectUnpack(args, v5Prime, v3Prime, dnCounter);
    }
    else if (args.specifiedBam == true){

        bamUnpack(args, v5Prime, v3Prime, dnCounter);        
    }

    //left forks
    std::vector<double> lRunOffs; 
    bool runLeft = false;
    bool runRight = false;

    if (args.specifiedLeft == true) {

        runLeft = true;

        simulation(v5Prime, v3Prime, lForkLength, lStallScore, args.specifiedDetect, args.specifiedBam, runLeft = true, runRight = false, lRunOffs);        
    }

    //right forks
    std::vector<double> rRunOffs;

    if (args.specifiedRight == true) {

        runRight = true;

        simulation(v5Prime, v3Prime, rForkLength, rStallScore, args.specifiedDetect, args.specifiedBam, runLeft = false, runRight = true, rRunOffs);
    }

    // Combine Simulation Run Offs
    std::vector<double> totalSimRunOffs;

    if (args.specifiedLeft == true && args.specifiedRight == true) {

        totalSimRunOffs.insert(totalSimRunOffs.end(), lRunOffs.begin(), lRunOffs.end());
        totalSimRunOffs.insert(totalSimRunOffs.end(), rRunOffs.begin(), rRunOffs.end());
    }
    else if (args.specifiedLeft == true) {

        totalSimRunOffs = lRunOffs;
    }
    else if (args.specifiedRight == true) {

        totalSimRunOffs = rRunOffs;
    }

    // Left Fork Observations
    std::vector<double> lObsRunOffs;
    if (args.specifiedLeft == true) {
        
        observation(lStallScore, lObsRunOffs);
    }

    // Right Fork Observations
    std::vector<double> rObsRunOffs;
    if (args.specifiedRight == true) {
        
        observation(rStallScore, rObsRunOffs);
    }

    // Combine Run Offs
    std::vector<double> totalObsRunOffs;
    if (args.specifiedLeft == true && args.specifiedRight == true) {

        totalObsRunOffs.insert(totalObsRunOffs.end(), lObsRunOffs.begin(), lObsRunOffs.end());
        totalObsRunOffs.insert(totalObsRunOffs.end(), rObsRunOffs.begin(), rObsRunOffs.end());
    }
    else if (args.specifiedLeft == true) {

        totalObsRunOffs = lObsRunOffs;
    }
    else if (args.specifiedRight == true) {

        totalObsRunOffs = rObsRunOffs;
    }

    // Calculate Simulation mean and standard deviation
    double simMean = vectorMean(totalSimRunOffs);
    double simStdDev = vectorStdv(totalSimRunOffs, simMean);

    // Calculate Observation mean abd standard deviation  
    double obsMean = vectorMean(totalObsRunOffs);
    double obsStdDev = vectorStdv(totalObsRunOffs, obsMean);

    // difference between simulated and observed
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

    // Output summary statistics
    std::cout << "Simulation Mean: " << simMean << "\n";
    std::cout << "Simulation Standard Deviation: " << simStdDev << "\n";
    std::cout << "Observation Mean: " << obsMean << "\n";
    std::cout << "Observation Standard Deviation: " << obsStdDev << "\n";
    std::cout << "#DifferenceMean " << difMean << "\n";
    std::cout << "#DifferenceStdv " << difStdDev << "\n";
    std::cout << "#95ConfidenceInterval " << leftTail << " " << rightTail << "\n";

    // Output data
    if (args.specifiedOutput == true) {

        std::ofstream outFile(args.output);
        
        outFile << "#DetectFile" << args.output << "\n";
        outFile << "#Software " << std::string(getExePath()) << "\n";
        outFile << "#Version " << std::string(VERSION) << "\n";
        outFile << "#Commit " << std::string(getGitCommit()) << "\n";
        outFile << "#n " << lnCounter + rnCounter << "\n";
        outFile << "#ExpectedMean " << simMean << "\n";
        outFile << "#ExpectedStdv " << simStdDev << "\n";
        outFile << "#ObservedMean " << obsMean << "\n";
        outFile << "#ObservedStdv " << obsStdDev << "\n";
        outFile << "#DifferenceMean " << difMean << "\n";
        outFile << "#DifferenceStdv " << difStdDev << "\n";
        outFile << "#95ConfidenceInterval " << leftTail << " " << rightTail << "\n";

        outFile << "\nsimulation:\n";

        for (const auto& val : totalSimRunOffs) {

            outFile << val << "\n";
        }

        outFile << "\nobservation:\n";

        for (const auto& val : totalObsRunOffs) {

            outFile << val << "\n";
        }
        outFile.close();
    }

    return 0;
    
}
