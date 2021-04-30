//----------------------------------------------------------
// Copyright 2019-2020 University of Oxford
// Written by Michael A. Boemo (mb915@cam.ac.uk)
// This software is licensed under GPL-3.0.  You should have
// received a copy of the license with this software.  If
// not, please Email the author.
//----------------------------------------------------------

#define TEST_CLUSTERING 0
//#define TEST_COOLDOWN 1

#include <fstream>
#include "regions.h"
#include "common.h"
#include "data_IO.h"
#include "error_handling.h"
#include "trainGMM.h"
#include <cmath>
#include <math.h>
#include <algorithm>
#define _USE_MATH_DEFINES


static const char *help=
"regions: DNAscent executable that finds regions of analogue incorporation from the output of DNAscent detect.\n"
"To run DNAscent regions, do:\n"
"   DNAscent regions -d /path/to/output.detect -o /path/to/output.regions\n"
"Required arguments are:\n"
"  -d,--detect               path to output file from DNAscent detect,\n"
"  -o,--output               path to output directory for bedgraph files.\n"
"Optional arguments (if used with HMM-based detect) are:\n"
"     --threshold            probability above which a BrdU call is considered positive (default: 0.8),\n"
"  -c,--cooldown             minimum gap between positive analogue calls (default: 4),\n"
"  -r,--resolution           number of thymidines in a region (default is 100 bp),\n"
"  -p,--probability          override probability that a thymidine 6mer contains a BrdU (default: automatically calculated),\n"
"  -z,--zScore               override zScore threshold for BrdU call (default: automatically calculated).\n"
"Written by Michael Boemo, Department of Pathology, University of Cambridge.\n"
"Please submit bug reports to GitHub Issues (https://github.com/MBoemo/DNAscent/issues).";

 struct Arguments {

	std::string detectFilename;
	double probability, threshold, likelihood;
	bool overrideProb,overrideZ,overrideResolution;
	unsigned int resolution;
	int cooldown;
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

 	//defaults - we'll override these if the option was specified by the user
	args.resolution = 100;
	args.threshold = 0;
	args.overrideProb = false;
	args.overrideZ = false;
	args.likelihood = 0.8;
	args.cooldown = 4;
	args.overrideResolution = false;

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
			args.overrideProb = true;
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
			args.overrideZ = true;
			i+=2;
		}
		else if ( flag == "--replication" ){

			args.callReplication = true;
			i+=1;
		}
		else if ( flag == "-c" or flag == "--cooldown" ){
 			std::string strArg( argv[ i + 1 ] );
			args.cooldown = std::stoi( strArg.c_str() );
			i+=2;
		}
		else if ( flag == "--threshold" ){
 			std::string strArg( argv[ i + 1 ] );
			args.likelihood = std::stof( strArg.c_str() );
			i+=2;
		}
		else if ( flag == "-r" or flag == "--resolution" ){
 			std::string strArg( argv[ i + 1 ] );
			args.resolution = std::stoi( strArg.c_str() );
			args.overrideResolution = true;
			i+=2;
		}
		else throw InvalidOption( flag );
	}
	if (args.outputFilename == args.detectFilename) throw OverwriteFailure();

	return args;
}


struct region{

	std::string call="";
	int start, end;
	double score;
	std::string forkDir="";
};


std::pair< double, double > twoMeans( std::vector< double > &observations ){

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

		//guard against either C1 or C2 being empty
		if (C1_points_new.size() == 0){

			C2_new = vectorMean(C2_points_new);
			return std::make_pair(0,C2_new);
		}
		else if (C2_points_new.size() == 0){

			C1_new = vectorMean(C1_points_new);
			return std::make_pair(C1_new,0);
		}

		C1_new = vectorMean(C1_points_new);
		C2_new = vectorMean(C2_points_new);

		C1_points_old = C1_points_new;
		C2_points_old = C2_points_new;

		iter++;
	}while (iter < maxIter and (std::abs(C1_old - C1_new)>tol or std::abs(C2_old - C2_new)>tol));

#if TEST_CLUSTERING
	std::cerr << ">" << C1_new << std::endl;
	for (auto c = C1_points_old.begin(); c < C1_points_old.end(); c++) std::cerr << *c << std::endl;
	std::cerr << ">" << C2_new << std::endl;
	for (auto c = C2_points_old.begin(); c < C2_points_old.end(); c++) std::cerr << *c << std::endl;
#endif 

	return std::make_pair(C1_new,C2_new);
}


int parseDetectLine_CNN(std::string line,
		        double callThreshold,
			unsigned int cooldownThreshold,
			 unsigned int &attemptCooldown,
			 unsigned int &callCooldown,
			 unsigned int &calls,
			 unsigned int &attempts){

	std::string column;
	std::stringstream ssLine(line);
	int position = -1, cIndex = 0;
	AnalogueScore B;
	while ( std::getline( ssLine, column, '\t' ) ){

		if ( cIndex == 0 ){

			position = std::stoi(column);
		}
		else if ( cIndex == 1 ){

			B.set(std::stof(column));
		}
		cIndex++;
	}
	assert(position != -1);


	if ( B.get() > callThreshold and position - callCooldown >= cooldownThreshold ){
		attemptCooldown = position;
		callCooldown = position;
		calls++;
		attempts++;
	}
	else if (position - attemptCooldown >= cooldownThreshold){
		attempts++;
		attemptCooldown = position;
	}
	return position;
}


std::string getStrand(std::string line){

	std::stringstream ssLine(line);
	int cIndex = 0;
	std::string strand = "";
	std::string column;
	int countCol = std::count(line.begin(), line.end(), ' ');
	assert(countCol == 4);
	while ( std::getline( ssLine, column, ' ' ) ){

		if ( cIndex == 4 ){

			strand = column;
		}
		cIndex++;
	}
	assert(strand != "");

	return strand;
}


void regionsCNN(Arguments args){

	//get a read count
	int readCount = 0;
	std::string line;
	std::ifstream inFile( args.detectFilename );
	if ( not inFile.is_open() ) throw IOerror( args.detectFilename );
	while( std::getline( inFile, line ) ){

		if ( line.substr(0,1) == ">" ) readCount++;
	}	
	progressBar pb(readCount,false);
	inFile.close();

	//estimate the fraction of BrdU incorporation
	double p;
	std::string header;
	unsigned int calls = 0, attempts = 0, gap = 0;

	int startingPos = -1;
	int progress = 0;
	unsigned int callCooldown = 0;
	unsigned int attemptCooldown = 0;
	std::string strand;

	if ( not args.overrideProb ){

	 	inFile.open( args.detectFilename );
		if ( not inFile.is_open() ) throw IOerror( args.detectFilename );

		std::cout << "Estimating analogue incorporation..." << std::endl;

		std::vector< double > callFractions;
		while( std::getline( inFile, line ) ){

			if (line.substr(0,1) == "#") continue; //ignore header
			if ( line.substr(0,1) == ">" ){

				strand = getStrand(line);
				progress++;
				pb.displayProgress( progress, 0, 0 );
				callCooldown = 0;
				attemptCooldown = 0;
				calls = 0, attempts = 0, gap = 0, startingPos = -1;
				continue;
			}

			int position = parseDetectLine_CNN(line, args.likelihood, args.cooldown, attemptCooldown, callCooldown, calls, attempts);

			if (position == -1) continue;

			if ( startingPos == -1 ) startingPos = position;
			gap = position - startingPos;

			if ( gap > args.resolution and attempts >= args.resolution / 10 ){

				double frac = (double) calls / (double) attempts;
				callFractions.push_back( frac );
				calls = 0, attempts = 0, gap = 0, startingPos = -1;
			}
		}

		std::cout << std::endl << "Done." << std::endl;
		double k1,k2;
		std::tie(k1,k2) = twoMeans( callFractions );
		p = std::max(k1,k2);

#if !TEST_CLUSTERING
		std::cerr << "Estimated fraction of analogue substitution in analogue-positive regions: " << p << std::endl;
#endif
		inFile.close();
	}
	else p = args.probability;

	if ( not args.overrideZ ){

		//estimate appropriate z-score threshold
		std::cout << "Setting Z-score threshold..." << std::endl;
		inFile.open( args.detectFilename );
		if ( not inFile.is_open() ) throw IOerror( args.detectFilename );
		progressBar pb_z(readCount,false);
		calls = 0; attempts = 0; gap = 0;
		startingPos = -1;
		progress = 0;
		callCooldown = 0; attemptCooldown = 0;
		std::vector<double> allZScores;
		while( std::getline( inFile, line ) ){

			if (line.substr(0,1) == "#") continue; //ignore header
			if ( line.substr(0,1) == ">" ){

				strand = getStrand(line);
				progress++;
				pb_z.displayProgress( progress, 0, 0 );
				calls = 0, attempts = 0, gap = 0, startingPos = -1;
				callCooldown = 0;
				attemptCooldown = 0;
			}
			else{

				int position = parseDetectLine_CNN(line, args.likelihood, args.cooldown, attemptCooldown, callCooldown, calls, attempts);

				if (position == -1) continue;
				if ( startingPos == -1 ) startingPos = position;
				gap = position - startingPos;

				if ( gap > args.resolution and attempts >= args.resolution / 10 ){

					double score = (calls - attempts * p) / sqrt( attempts * p * ( 1 - p) );
					allZScores.push_back(score); 
					calls = 0, attempts = 0, gap = 0, startingPos = -1;
				}
			}
		}
		inFile.close();
		std::cout << "Done." << std::endl;

		std::vector< double > fitParams = gaussianMixtureEM(0.5, -7.0, 1.0, 0., 1.0, allZScores, 0.01, 100 );
		double thym_mu, thym_mix, thym_sigma, brdu_mu, brdu_mix, brdu_sigma;
		if (fitParams[1] < fitParams[2]){

			thym_mix = fitParams[0];
			thym_mu = fitParams[1];
			thym_sigma = fitParams[2];
			brdu_mix = fitParams[3];
			brdu_mu = fitParams[4];
			brdu_sigma = fitParams[5];
		}
		else{

			thym_mix = fitParams[3];
			thym_mu = fitParams[4];
			thym_sigma = fitParams[5];
			brdu_mix = fitParams[0];
			brdu_mu = fitParams[1];
			brdu_sigma = fitParams[2];
		}
#if !TEST_CLUSTERING
		std::cerr << "Estimated fraction of thymidine regions: " << thym_mix << std::endl;
		std::cerr << "Estimated fraction of BrdU regions: " << brdu_mix << std::endl;
		std::cerr << "Thymidine Z-score mean, stdv: " << thym_mu << " " << thym_sigma << std::endl;
		std::cerr << "BrdU Z-score mean, stdv: " << brdu_mu << " " << brdu_sigma << std::endl;
#endif

#if !TEST_CLUSTERING
		std::cerr << "Set Z-score threshold: " << brdu_mu << std::endl;
#endif
		args.threshold = brdu_mu;
	}
	
	//call regions
 	inFile.open( args.detectFilename );
	if ( not inFile.is_open() ) throw IOerror( args.detectFilename );
 	std::ofstream outFile( args.outputFilename );
	if ( not outFile.is_open() ) throw IOerror( args.outputFilename );

	//write the regions header
	outFile <<  writeRegionsHeader(args.detectFilename, args.likelihood, true, args.cooldown, args.resolution, p, args.threshold);

	std::ofstream repFile;
	if ( args.callReplication ){

		repFile.open("calledOrigins.dnascent");
	}
	std::cout << "Calling regions..." << std::endl;
	std::vector< region > buffer;
	calls = 0; attempts = 0; gap = 0;
	startingPos = -1;
	progress = 0;
	callCooldown = 0; attemptCooldown = 0;
	bool first = true;
	while( std::getline( inFile, line ) ){

		if (line.substr(0,1) == "#") continue; //ignore header
		if ( line.substr(0,1) == ">" ){

			strand = getStrand(line);
			progress++;
			pb.displayProgress( progress, 0, 0 );

			if (not first){

				outFile << header << std::endl;
				
				for ( auto r = buffer.begin(); r < buffer.end(); r++ ){

					outFile << r -> start << "\t" << r -> end << "\t" << r -> score << "\t" << r -> call << "\t" << r -> forkDir <<  std::endl;
				}
			}
			header = line;
			buffer.clear();
			calls = 0, attempts = 0, gap = 0, startingPos = -1;
			callCooldown = 0;
			attemptCooldown = 0;
			first = false;
			
		}
		else{

			int position = parseDetectLine_CNN(line, args.likelihood, args.cooldown, attemptCooldown, callCooldown, calls, attempts);

			if (position == -1) continue;

			if ( startingPos == -1 ) startingPos = position;
			gap = position - startingPos;

			if ( gap > args.resolution and attempts >= args.resolution / 10 ){

				region r;

				r.score = (calls - attempts * p) / sqrt( attempts * p * ( 1 - p) );

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
				
		for ( auto r = buffer.begin(); r < buffer.end(); r++ ){

			outFile << r -> start << "\t" << r -> end << "\t" << r -> score << "\t" << r -> call << "\t" << r -> forkDir <<  std::endl;
		}
	}

	if ( repFile.is_open() ) repFile.close();
	inFile.close();
	outFile.close();
	std::cout << std::endl << "Done." << std::endl;

	if (p < 0.1){
		std::cerr << "WARNING: Analogue incorporation is estimated to be low: " << p << std::endl;
		std::cerr << "   Samples may not have analogue in them; DNAscent regions assumes there are both analogue-positive and analogue-negative regions in the sample." << std::endl;
		std::cerr << "   The DNAscent regions results may be unreliable. See https://dnascent.readthedocs.io/en/latest/regions.html for details." << std::endl;
	}
	else if (p > 0.7){
		std::cerr << "WARNING: Analogue incorporation is estimated to be high: " << p << std::endl;
		std::cerr << "   Samples may be saturated; DNAscent regions assumes there are both analogue-positive and analogue-negative regions in the sample." << std::endl;
		std::cerr << "   The DNAscent regions results may be unreliable. See https://dnascent.readthedocs.io/en/latest/regions.html for details." << std::endl;
	}
}


int regions_main( int argc, char** argv ){

	Arguments args = parseRegionsArguments( argc, argv );

	regionsCNN(args);

	return 0;
}
