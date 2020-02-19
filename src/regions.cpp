//----------------------------------------------------------
// Copyright 2017 University of Oxford
// Written by Michael A. Boemo (mb915@cam.ac.uk)
// This software is licensed under GPL-2.0.  You should have
// received a copy of the license with this software.  If
// not, please Email the author.
//----------------------------------------------------------

#define TEST_CLUSTERING 0
//#define TEST_COOLDOWN 1

#include <fstream>
#include "regions.h"
#include "data_IO.h"
#include "error_handling.h"
#include "train.h"
#include <cmath>
#include <math.h>
#include <algorithm>
#define _USE_MATH_DEFINES


static const char *help=
"regions: DNAscent executable that finds regions of analogue incorporation from the output of DNAscent detect.\n"
"To run DNAscent regions, do:\n"
"  ./DNAscent regions -d /path/to/output.detect -o /path/to/output.regions\n"
"Required arguments are:\n"
"  -d,--detect               path to output file from DNAscent detect,\n"
"  -o,--output               path to output directory for bedgraph files.\n"
"Optional arguments are:\n"
"     --replication          detect fork direction and call origin firing (default: off),\n"
"  -l,--likelihood           log-likelihood threshold for a positive analogue call (default: 1.25),\n"
"  -c,--cooldown             minimum gap between positive analogue calls (default: 4),\n"
"  -r,--resolution           minimum length of regions (default is 2kb).\n"
"  -p,--probability          override probability that a thymidine 6mer contains a BrdU (default: automatically calculated),\n"
"  -z,--zScore               override zScore threshold for BrdU call (default: automatically calculated).\n"
"Written by Michael Boemo, Department of Pathology, University of Cambridge.\n"
"Please submit bug reports to GitHub Issues (https://github.com/MBoemo/DNAscent/issues).";

 struct Arguments {

	std::string detectFilename;
	double probability, threshold, likelihood;
	bool overrideProb,overrideZ;
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
	args.resolution = 2000;
	args.threshold = 0;
	args.callReplication = false;
	args.overrideProb = false;
	args.overrideZ = false;
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
		else if ( flag == "-l" or flag == "--likelihood" ){
 			std::string strArg( argv[ i + 1 ] );
			args.likelihood = std::stof( strArg.c_str() );
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

	std::string call="";
	int start, end;
	double score;
	std::string forkDir="";
};


double vectorMean( std::vector< double > &obs ){

	double total = 0.0;
	for ( size_t i = 0; i < obs.size(); i++ ) total += obs[i];
	return total / (double) obs.size();
}


std::pair< double, double > twoMeans( std::vector< double > &observations ){

	double C1_old = 0.01;
	double C2_old = 0.3;
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

#if TEST_CLUSTERING
	std::cerr << ">" << C1_new << std::endl;
	for (auto c = C1_points_old.begin(); c < C1_points_old.end(); c++) std::cerr << *c << std::endl;
	std::cerr << ">" << C2_new << std::endl;
	for (auto c = C2_points_old.begin(); c < C2_points_old.end(); c++) std::cerr << *c << std::endl;
#endif 

	return std::make_pair(C1_new,C2_new);
}


void callOrigins( std::vector< region > &regions, double threshold, std::string &header, std::ofstream &repFile ){

	//10kb moving average filter
	std::vector< double > out;
	for ( unsigned int i = 2; i < regions.size() - 2; i++ ){

		out.push_back( ( regions[i-2].score + regions[i-1].score + regions[i].score + regions[i+1].score + regions[i+2].score ) / 5.0 );
	}
	for ( unsigned int i = 0; i < out.size(); i++ ){

		regions[i+2].score = out[i];
	}

	//secondary 6kb moving average filter
	out.clear();
	for ( unsigned int i = 1; i < regions.size() - 1; i++ ){

		out.push_back( ( regions[i-1].score + regions[i].score + regions[i+1].score ) / 3.0 );
	}
	for ( unsigned int i = 0; i < out.size(); i++ ){

		regions[i+1].score = out[i];
		
		//re-assess the call with the smoothed scores
		if ( regions[i+1].score >= 0 ) regions[i+1].call = "BrdU";
		else regions[i+1].call = "Thym";
	}
	out.clear();

	std::vector< double > derivatives;
	
	//get the central derivative for all the middle regions
	for ( unsigned int i = 1; i < regions.size() - 1; i++ ){

		double former = regions[i-1].score;
		double next = regions[i+1].score;

		derivatives.push_back( ( next - former ) / 2.0 );
	}

	//error correct for derivatives that are sandwiched between two of the opposite sign
	for ( unsigned int i = 2; i < derivatives.size() - 2; i++ ){

		if ( derivatives[i-1] < 0 and derivatives[i+1] < 0 and derivatives[i] > 0 ) derivatives[i] = ( derivatives[i-1] + derivatives[i+1] ) / 2.0;
		if ( derivatives[i-1] > 0 and derivatives[i+1] > 0 and derivatives[i] < 0 ) derivatives[i] = ( derivatives[i-1] + derivatives[i+1] ) / 2.0;
	}

	//smooth derivative calls
	for ( unsigned int i = 1; i < derivatives.size() - 1; i++ ){

		out.push_back( ( derivatives[i-1] + derivatives[i] + derivatives[i+1] ) / 3.0 );
	}
	for ( unsigned int i = 0; i < out.size(); i++ ){

		derivatives[i+1] = out[i];
	}
	
	//error correct on BrdU calls for the middle
	for ( unsigned int i = 1; i < regions.size() - 1; i++ ){

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
	for ( unsigned int i = 1; i < regions.size() - 1; i++ ){

		if ( regions[i].call == "BrdU" and derivatives[i] < 0 ) regions[i].forkDir = "+";
		if ( regions[i].call == "BrdU" and derivatives[i] > 0 ) regions[i].forkDir = "-";
	}

	//call origin positions
	std::vector< int > oriPosForRead;
	for ( unsigned int i = 1; i < regions.size()-2; i++ ){

		if ( regions[i-1].forkDir == "-" and regions[i].forkDir == "-" and regions[i+1].forkDir == "+" and regions[i+2].forkDir == "+" ){

			if ( regions[i].score > regions[i+1].score ) oriPosForRead.push_back( (regions[i].start + regions[i].end) / 2 );
			else oriPosForRead.push_back( (regions[i+1].start + regions[i+1].end) / 2 );
		}
	}

	//fix the ends
	if (regions.size() > 3){

		if (not regions[1].forkDir.empty() and regions[0].call == "BrdU") regions[0].forkDir = regions[1].forkDir;
		if (not regions[regions.size()-2].forkDir.empty() and regions[regions.size()-1].call == "BrdU") regions[regions.size()-1].forkDir = regions[regions.size()-2].forkDir;
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
	int callCooldown = 0;
	int attemptCooldown = 0;

	if ( not args.overrideProb ){

	 	inFile.open( args.detectFilename );
		if ( not inFile.is_open() ) throw IOerror( args.detectFilename );

		std::cout << "Estimating analogue incorporation..." << std::endl;

		std::vector< double > callFractions;
		while( std::getline( inFile, line ) ){

			if ( line.substr(0,1) == ">" ){

				progress++;
				pb.displayProgress( progress, 0, 0 );
				callCooldown = 0;
				attemptCooldown = 0;
				continue;
			}

			std::string column;
			std::stringstream ssLine(line);
			int position = -1, cIndex = 0;
			AnalogueScore B, BM;
			int countCol = std::count(line.begin(), line.end(), '\t');
			while ( std::getline( ssLine, column, '\t' ) ){

				if ( cIndex == 0 ){

					position = std::stoi(column);
				}
				else if ( cIndex == 1 ){

					B.set(std::stof(column));
				}
				else if ( cIndex == 2 and countCol > 3 ){ //methyl-aware detect file

					BM.set(std::stof(column));
				}
				cIndex++;
			}
			assert(position != -1);
			if ( countCol > 3){

					if ( B.get() > args.likelihood and BM.get() > args.likelihood ){
						calls++;
						attempts++;
					}
					else if ( B.get() < args.likelihood and BM.get() > 0.0 ) attempts++;
					// if B < args.likelihood and BM < 0, then it's corrupted by methylation so don't count as attempt
			}
			else{
					//testing
					//std::cout << B.get() << " " << position << " " << callCooldown << std::endl;
					if ( B.get() > args.likelihood and position - callCooldown >= args.cooldown ){
						attemptCooldown = position;
						callCooldown = position;
						calls++;
						attempts++;
					}
					else if (position - attemptCooldown >= args.cooldown){
						attempts++;
						attemptCooldown = position;
					}
			}

			if ( startingPos == -1 ) startingPos = position;
			gap = position - startingPos;

			if ( gap > args.resolution and attempts >= args.resolution / 30 ){

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

			if ( line.substr(0,1) == ">" ){

				progress++;
				pb_z.displayProgress( progress, 0, 0 );
				calls = 0, attempts = 0, gap = 0, startingPos = -1;
				callCooldown = 0;
				attemptCooldown = 0;
			}
			else{

				std::string column;
				std::stringstream ssLine(line);
				int position = -1, cIndex = 0;
				AnalogueScore B, BM;
				int countCol = std::count(line.begin(), line.end(), '\t');
				while ( std::getline( ssLine, column, '\t' ) ){

					if ( cIndex == 0 ){

						position = std::stoi(column);
					}
					else if ( cIndex == 1 ){

						B.set(std::stof(column));
					}
					else if ( cIndex == 2 and countCol > 3 ){ //methyl-aware detect file

						BM.set(std::stof(column));
					}
					cIndex++;
				}
				assert(position != -1);

				if ( countCol > 3){

						if ( B.get() > args.likelihood and BM.get() > args.likelihood ){
							calls++;
							attempts++;
						}
						else if ( B.get() < args.likelihood and BM.get() > 0 ) attempts++;
						//not strong BrdU but more BrdU than methyl counts as an attempt
				}
				else{

					if ( B.get() > args.likelihood and position - callCooldown >= args.cooldown){
						callCooldown = position;
						attemptCooldown = position;
						calls++;
						attempts++;
					}
					else if (position - attemptCooldown >= args.cooldown){
						attempts++;
						attemptCooldown = position;
					}
				}

				if ( startingPos == -1 ) startingPos = position;
				gap = position - startingPos;

				if ( gap > args.resolution and attempts >= args.resolution / 30 ){

					double score = (calls - attempts * p) / sqrt( attempts * p * ( 1 - p) );
					allZScores.push_back(score);
					calls = 0, attempts = 0, gap = 0, startingPos = -1;
				}
			}
		}
		inFile.close();
		std::cout << "Done." << std::endl;

		std::vector< double > fitParams = gaussianMixtureEM(-3.0, 3.0, 0, 3.0, allZScores, 0.01, 100 );
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
		if (2*thym_sigma < (brdu_mu - thym_mu)/2.0){

#if !TEST_CLUSTERING
			std::cerr << "Set Z-score threshold: " << thym_mu+(brdu_mu - thym_mu)/2.0 << std::endl;
#endif
			args.threshold = thym_mu+(brdu_mu - thym_mu)/2.0;
		}
		else{

#if !TEST_CLUSTERING
			std::cerr << "Set Z-score threshold: " << thym_mu+2*thym_sigma << std::endl;
#endif
			args.threshold = thym_mu+2*thym_sigma;
		}
	}
	
	//call regions
 	inFile.open( args.detectFilename );
	if ( not inFile.is_open() ) throw IOerror( args.detectFilename );
 	std::ofstream outFile( args.outputFilename );
	if ( not outFile.is_open() ) throw IOerror( args.outputFilename );

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
	while( std::getline( inFile, line ) ){

		if ( line.substr(0,1) == ">" ){

			progress++;
			pb.displayProgress( progress, 0, 0 );

			if ( buffer.size() > 5 and args.callReplication ){

				outFile << header << std::endl;
			
				callOrigins(buffer,args.threshold, header, repFile);
				
				for ( auto r = buffer.begin(); r < buffer.end(); r++ ){

					outFile << r -> start << "\t" << r -> end << "\t" << r -> score << "\t" << r -> call << "\t" << r -> forkDir <<  std::endl;
				}
			}
			else if (not args.callReplication){

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
			
		}
		else{

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

			if ( startingPos == -1 ) startingPos = position;
			gap = position - startingPos;

			if ( B.get() > args.likelihood and position - callCooldown >= args.cooldown){
#if TEST_COOLDOWN
std::cout << ">>>>>Call: " << line << std::endl;
std::cout << "     Call Cooldown: " << callCooldown << std::endl;
std::cout << "     Attempt Cooldown: " << attemptCooldown << std::endl;
#endif
				callCooldown = position;
				attemptCooldown = position;
				calls++;
				attempts++;
			}
			else if (position - attemptCooldown >= args.cooldown){
#if TEST_COOLDOWN
std::cout << ">>>>>Attempt: " << line << std::endl;
std::cout << "     Call Cooldown: " << callCooldown << std::endl;
std::cout << "     Attempt Cooldown: " << attemptCooldown << std::endl;
#endif
				attempts++;
				attemptCooldown = position;
			}
#if TEST_COOLDOWN
			else{
			std::cout << ">>>>>Ignored: " << line << std::endl;
			std::cout << "     Call Cooldown: " << callCooldown << std::endl;
			std::cout << "     Attempt Cooldown: " << attemptCooldown << std::endl;
			}
#endif
			if ( gap > args.resolution and attempts >= args.resolution / 30 ){

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
			
		if (args.callReplication) callOrigins(buffer,args.threshold, header, repFile);
				
		for ( auto r = buffer.begin(); r < buffer.end(); r++ ){

			outFile << r -> start << "\t" << r -> end << "\t" << r -> score << "\t" << r -> call << "\t" << r -> forkDir <<  std::endl;
		}
	}

	if ( repFile.is_open() ) repFile.close();
	inFile.close();
	outFile.close();
	std::cout << std::endl << "Done." << std::endl;

	return 0;
}
