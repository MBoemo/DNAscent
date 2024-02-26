//----------------------------------------------------------
// Copyright 2019 University of Oxford
// This software is licensed under GPL-3.0.  You should have
// received a copy of the license with this software.  If
// not, please Email the author.
//----------------------------------------------------------

#include <exception>
#include <algorithm>
#include <math.h>
#include <vector>
#include <iostream>
#include <fstream>
#include <cmath>
#include "common.h"
#include "data_IO.h"
#include "error_handling.h"
#include "event_handling.h"
#include "probability.h"
#include "trainGMM.h"
#include "config.h"

static const char *help=
"trainGMM: DNAscent executable that determines the mean and standard deviation of a base analogue's current.\n"
"Note: This executable is geared towards developers and advanced users.\n"
"To run DNAscent trainGMM, do:\n"
"   DNAscent trainGMM -d /path/to/DNAscent.align -o output.model\n"
"Required arguments are:\n"
"  -d,--trainingData         path to training data from DNAscent align,\n"
"  -o,--output               path to the output trained pore model file.\n"
"Optional arguments are:\n"
"  -pi,                      mixing parameter for BrdU (default is 0.5),\n"
"  -m,--max-reads            maximum number of reads to consider (default is 100000),\n"
"  -e,--max-events           maximum number of events per 6mer to consider (default is 10000),\n"
"  -t,--threads              number of threads (default is 1 thread).\n"
"DNAscent is under active development by the Boemo Group, Department of Pathology, University of Cambridge (https://www.boemogroup.org/).\n"
"Please submit bug reports to GitHub Issues (https://github.com/MBoemo/DNAscent/issues).";

struct Arguments {

	std::string trainingOutputFilename;
	bool logFile;
	std::string logFilename;
	int threads;
	std::string eventalignFilename;
	unsigned int maxReads;
	bool capReads;
	unsigned int maxEvents;
	float pi;
};

Arguments parseTrainingArguments( int argc, char** argv ){

	if( argc < 2 ){

		std::cout << "Exiting with error.  Insufficient arguments passed to DNAscent train." << std::endl << help << std::endl;
		exit(EXIT_FAILURE);
	}

	if ( std::string( argv[ 1 ] ) == "-h" or std::string( argv[ 1 ] ) == "--help" ){

		std::cout << help << std::endl;
		exit(EXIT_SUCCESS);
	}

	Arguments trainArgs;

	//defaults - we'll override these if the option was specified by the user
	trainArgs.threads = 1;
	trainArgs.maxReads = 100000;
	trainArgs.capReads = false;
	trainArgs.maxEvents = 10000;
	trainArgs.pi = 0.5;

	/*parse the command line arguments */
	for ( int i = 1; i < argc; ){

		std::string flag( argv[ i ] );

		if ( flag == "-d" or flag == "--trainingData" ){

			std::string strArg( argv[ i + 1 ] );
			trainArgs.eventalignFilename = strArg;
			i+=2;
		}
		else if ( flag == "-o" or flag == "--output" ){

			std::string strArg( argv[ i + 1 ] );
			trainArgs.trainingOutputFilename = strArg;
			i+=2;
		}
		else if ( flag == "-m" or flag == "--max-reads" ){

			std::string strArg( argv[ i + 1 ] );
			trainArgs.maxReads = std::stoi(strArg.c_str());
			trainArgs.capReads = true;
			i+=2;
		}
		else if ( flag == "-e" or flag == "--max-events" ){

			std::string strArg( argv[ i + 1 ] );
			trainArgs.maxEvents = std::stoi(strArg.c_str());
			i+=2;
		}
		else if ( flag == "-pi" ){

			std::string strArg( argv[ i + 1 ] );
			trainArgs.pi = std::stof(strArg.c_str());
			i+=2;
		}
		else if ( flag == "-o" or flag == "--output" ){

			std::string strArg( argv[ i + 1 ] );
			trainArgs.trainingOutputFilename = strArg;
			i+=2;
		}
		else if ( flag == "-t" or flag == "--threads" ){

			std::string strArg( argv[ i + 1 ] );
			trainArgs.threads = std::stoi( strArg.c_str() );
			i+=2;
		}
		else throw InvalidOption( flag );
	}
	if (trainArgs.trainingOutputFilename == trainArgs.eventalignFilename) throw OverwriteFailure();

	return trainArgs;
}


std::vector<int> findNeighbours( std::vector<double> &events, double ev, double epsilon ){

	std::vector< int > neighbourIdx;
	for ( unsigned int i = 0; i < events.size(); i++ ){

		if (std::abs(ev - events[i]) <= epsilon) neighbourIdx.push_back(i);
	}
	return neighbourIdx;
}


std::map<int,int> DBSCAN( std::vector< double > &events, double epsilon, unsigned int minPoints ){

	//labels
	//-2 := undefined
	//-1 := noise
	// 0 <= cluster int

	//initialise labels
	std::map< int, int > index2label;
	for ( unsigned int i = 0; i < events.size(); i++ ) index2label[i] = -2;

	int C = 0; //cluster counter
	for ( unsigned int i = 0; i < events.size(); i++ ){

		if (index2label[i] != -2) continue;
		std::vector<int> neighbourIndices = findNeighbours( events, events[i], epsilon );
		if (neighbourIndices.size() < minPoints) {

			index2label[i] = -1; //label as noise
			continue;
		}

		C++; //increment the cluster
		index2label[i] = C;
		std::vector< int > seedSet = neighbourIndices;
		seedSet.erase(std::find(seedSet.begin(),seedSet.end(),i)); //seed set is the neighbours minus the event we're at
		for ( unsigned int j = 0; j < seedSet.size(); j++ ){

			if (index2label[seedSet[j]] == -1) index2label[seedSet[j]] = C;
			if (index2label[seedSet[j]] != -2 ) continue;
			index2label[seedSet[j]] = C;
			std::vector<int> neighbourIndicesInner = findNeighbours(events, events[seedSet[j]], epsilon);
			if (neighbourIndicesInner.size() >= minPoints){

				seedSet.insert(seedSet.end(),neighbourIndicesInner.begin(),neighbourIndicesInner.end());
			}
		}
	}
	return index2label;
}


std::vector< double > gaussianMixtureEM_PRIOR( double pi, double mu1, double sigma1, double mu2, double sigma2, std::vector< double > &data, double tolerance, int maxIter ){
//only trains N(mu2,sigma2) and leaves N(mu1,sigma1) alone

	std::vector< std::vector< double > > Z( 2, std::vector< double >( data.size(), 0 ) );

	double pi1 = 1-pi;
	double pi2 = pi;
	double total1, total2, logLikelihood_Old, logLikelihood_New;

	/*INITIALISATION - calculuate the log likelihood using the paramters initially passed */
	logLikelihood_Old = 0;
	double expll;
	for ( unsigned int i = 0; i < data.size(); i++ ){

		expll = pi1 * normalPDF(mu1,sigma1,data[i]) + pi2 * normalPDF(mu2,sigma2,data[i]);
		logLikelihood_Old += eln(expll);
	}
	double improvement = std::numeric_limits< double >::max();

	int iterations = 0;
	while ( improvement > tolerance ){

		/*EXPECTATION */
		for ( unsigned int i = 0; i < data.size(); i++ ){

			/*first normal distribution */
			Z[0][i] = pi1 * normalPDF(mu1,sigma1,data[i]) / (pi1 * normalPDF(mu1,sigma1,data[i]) + pi2 * normalPDF(mu2,sigma2,data[i]));

			/*second normal distribution */
			Z[1][i] = pi2 * normalPDF(mu2,sigma2,data[i]) / (pi1 * normalPDF(mu1,sigma1,data[i]) + pi2 * normalPDF(mu2,sigma2,data[i]));
		}

		/*MAXIMISATION */
		/*calculuate Nk's from Z */
		double Nk1 = 0;
		double Nk2 = 0;
		for ( unsigned int i = 0; i < data.size(); i++ ){

			Nk1 += Z[0][i];
			Nk2 += Z[1][i];
		}

		/*re-estimate pi */
		pi1 = Nk1 / (double) data.size();
		pi2 = Nk2 / (double) data.size();

		/*re-estimate mu */
		total1 = 0;
		total2 = 0;
		for ( unsigned int i = 0; i < data.size(); i++ ){

			total1 += Z[0][i]*data[i];
			total2 += Z[1][i]*data[i];
		}
		mu2 = total2 / Nk2;

		/*re-estimate sigma */
		total1 = 0;
		total2 = 0;
		for ( unsigned int i = 0; i < data.size(); i++ ){

			total1 += Z[0][i]*(data[i] - mu1)*(data[i] - mu1);
			total2 += Z[1][i]*(data[i] - mu2)*(data[i] - mu2);
		}
		sigma2 = sqrt( total2 / Nk2 );

		/*compute new likelihood */
		logLikelihood_New = 0;
		for ( unsigned int i = 0; i < data.size(); i++ ){

			expll = pi1 * normalPDF(mu1,sigma1,data[i]) + pi2 * normalPDF(mu2,sigma2,data[i]);
			logLikelihood_New += eln(expll);
		}

		improvement = logLikelihood_New - logLikelihood_Old;
		logLikelihood_Old = logLikelihood_New;
		iterations++;
		if (iterations > maxIter) break;
	}
	return { pi1, mu1, sigma1, pi2, mu2, sigma2 };
}


std::vector< double > gaussianMixtureEM( double pi, double mu1, double sigma1, double mu2, double sigma2, std::vector< double > &data, double tolerance, int maxIter ){

	std::vector< std::vector< double > > Z( 2, std::vector< double >( data.size(), 0 ) );

	double pi1 = 1-pi;
	double pi2 = pi;
	double total1, total2, logLikelihood_Old, logLikelihood_New;

	/*INITIALISATION - calculuate the log likelihood using the paramters initially passed */
	logLikelihood_Old = 0;
	double expll;
	for ( unsigned int i = 0; i < data.size(); i++ ){

		expll = pi1 * normalPDF(mu1,sigma1,data[i]) + pi2 * normalPDF(mu2,sigma2,data[i]);
		logLikelihood_Old += eln(expll);
	}
	double improvement = std::numeric_limits< double >::max();

	int iterations = 0;
	while ( improvement > tolerance ){

		/*EXPECTATION */
		for ( unsigned int i = 0; i < data.size(); i++ ){

			/*first normal distribution */
			Z[0][i] = pi1 * normalPDF(mu1,sigma1,data[i]) / (pi1 * normalPDF(mu1,sigma1,data[i]) + pi2 * normalPDF(mu2,sigma2,data[i]));

			/*second normal distribution */
			Z[1][i] = pi2 * normalPDF(mu2,sigma2,data[i]) / (pi1 * normalPDF(mu1,sigma1,data[i]) + pi2 * normalPDF(mu2,sigma2,data[i]));
		}

		/*MAXIMISATION */
		/*calculuate Nk's from Z */
		double Nk1 = 0;
		double Nk2 = 0;
		for ( unsigned int i = 0; i < data.size(); i++ ){

			Nk1 += Z[0][i];
			Nk2 += Z[1][i];
		}

		/*re-estimate pi */
		pi1 = Nk1 / (double) data.size();
		pi2 = Nk2 / (double) data.size();

		/*re-estimate mu */
		total1 = 0;
		total2 = 0;
		for ( unsigned int i = 0; i < data.size(); i++ ){

			total1 += Z[0][i]*data[i];
			total2 += Z[1][i]*data[i];
		}
		mu1 = total1 / Nk1;
		mu2 = total2 / Nk2;

		/*re-estimate sigma */
		total1 = 0;
		total2 = 0;
		for ( unsigned int i = 0; i < data.size(); i++ ){

			total1 += Z[0][i]*(data[i] - mu1)*(data[i] - mu1);
			total2 += Z[1][i]*(data[i] - mu2)*(data[i] - mu2);
		}
		sigma1 = sqrt( total1 / Nk1 );
		sigma2 = sqrt( total2 / Nk2 );

		/*compute new likelihood */
		logLikelihood_New = 0;
		for ( unsigned int i = 0; i < data.size(); i++ ){

			expll = pi1 * normalPDF(mu1,sigma1,data[i]) + pi2 * normalPDF(mu2,sigma2,data[i]);
			logLikelihood_New += eln(expll);
		}

		improvement = logLikelihood_New - logLikelihood_Old;
		logLikelihood_Old = logLikelihood_New;
		iterations++;
		if (iterations > maxIter) break;
	}
	return { pi1, mu1, sigma1, pi2, mu2, sigma2 };
}


void printAllKLengthRec(char set[], std::string prefix, int n, int k, std::vector<std::string> &out){

	//base case
	if (k == 0){
		out.push_back(prefix);
		return;
	}

	//recursion 
	for (int i = 0; i < n; i++){

		std::string newPrefix;
		newPrefix = prefix + set[i];
		printAllKLengthRec(set, newPrefix, n, k - 1, out);
	}
}


void printAllKLength(char set[], int k,int n, std::vector<std::string> &out){
    
	printAllKLengthRec(set, "", n, k, out);
}


int train_main( int argc, char** argv ){

	Arguments trainArgs = parseTrainingArguments( argc, argv );

	/*open output file */
	std::ofstream outFile( trainArgs.trainingOutputFilename );
	if ( not outFile.is_open() ) throw IOerror( trainArgs.trainingOutputFilename );

	std::string line;
	int prog, failed;

	/*fudge for openmp */
	char set1[] = {'A', 'T', 'G', 'C'};
	unsigned int k = Pore_Substrate_Config.kmer_len;
	std::vector<std::string> allKmers;
	printAllKLength(set1, k, 4, allKmers);

	std::map< int, std::string > intToKmer;
	std::map< std::string, int > kmerToInt;
	int index = 0;
	for ( unsigned int i = 0; i < allKmers.size(); i++ ){

		intToKmer[index] = allKmers[i];
		kmerToInt[allKmers[i]] = index;
		index++;
	}

	std::vector< std::vector< double > > importedEvents( pow(4,k) );

	//get a read count
	unsigned int readCount = 0;

	if (not trainArgs.capReads){
		std::ifstream readStream(trainArgs.eventalignFilename);	
		if ( not readStream.is_open() ) throw IOerror( trainArgs.eventalignFilename );
		while( std::getline( readStream, line ) ){
			if ( line.substr(0,1) == ">" ) readCount++;
		}
		readStream.close();	
	
	}
	else{
		readCount = trainArgs.maxReads;
	}
	
	progressBar pb_read(readCount,true);

	unsigned int readsRead = 0;
	std::ifstream eventFile(trainArgs.eventalignFilename);	
	if ( not eventFile.is_open() ) throw IOerror( trainArgs.eventalignFilename );
	while ( std::getline( eventFile, line) ){
	
		if (line.empty()) continue;

		if ( line.substr(0,1) == "#" ) continue;

		if ( line.substr(0,1) == ">" ){

			readsRead++;
			pb_read.displayProgress( readsRead, 0, 0 );
			continue;
		}

		std::istringstream ss( line );
		std::string kmer, entry;
		double eventMean;

		int col = 0;
		while ( std::getline( ss, entry, '\t' ) ){

			if ( col == 3 ){

				kmer = entry;
				break;
			}
			else if ( col == 2 ){

				eventMean = atof( entry.c_str() );
			}
			col++;
		}

		if ( importedEvents[kmer2index(kmer, k)].size() < trainArgs.maxEvents ){
			importedEvents[kmer2index(kmer, k)].push_back( eventMean );
		}
		if (readsRead > trainArgs.maxReads) break;
	}
	eventFile.close();
	std::cout << "ok." << std::endl;
	std::cout << "Fitting..." << std::endl;
	
	/*fit a mixture model to the events that aligned to each position in the reference */
	outFile << "6mer" << '\t' << "ONT_mean" << '\t' << "ONT_stdv" << '\t' << "pi_1" << '\t' << "mean_1" << '\t' << "stdv_1" << '\t' << "pi_2" << '\t' << "mean_2" << '\t' << "stdv_2" << '\t' << "imported_events" << '\t' << "filtered_events" << std::endl;
	progressBar pb_fit( importedEvents.size(),true );

	#pragma omp parallel for schedule(dynamic) shared(pb_fit, Pore_Substrate_Config, prog, failed, outFile, importedEvents, trainArgs) num_threads(trainArgs.threads)
	for ( unsigned int i = 0; i < importedEvents.size(); i++ ){

		/*don't train if we have less than 200 events for this 6mer */
		if ( importedEvents[i].size() < 200 ){
			prog++;
			continue;
		}

		//DBSCAN to eliminate alignment artefacts
		unsigned int minPoints = 0.025*importedEvents[i].size();
		std::map<int,int> labels = DBSCAN( importedEvents[i], 0.5, minPoints );

		std::vector< double > filteredEvents;
		for ( unsigned int j = 0; j < importedEvents[i].size(); j++ ){

			if (labels[j] != -1) filteredEvents.push_back(importedEvents[i][j]);
		}

		if ( filteredEvents.size() < 50 ){
			prog++;
			continue;
		}

		std::string kmer = intToKmer[i];
		double mu1, stdv1, mu2, stdv2;

		/*get the ONT distribution for the mixture */
		std::pair<double,double> meanStd = Pore_Substrate_Config.pore_model[kmer2index(kmer, k)];
		mu1 = meanStd.first;
		stdv1 = meanStd.second;

		/*make a second distribution that's similar to the ONT distribution */
		mu2 = meanStd.first;
		stdv2 = 2*meanStd.second;

		/*fit the model */
		std::vector< double > fitParameters;
		try{

			fitParameters = gaussianMixtureEM_PRIOR( trainArgs.pi, mu1, stdv1, mu2, stdv2, filteredEvents, 0.01, 100 );
		}
		catch ( NegativeLog &nl ){

			failed++;
			prog++;
			continue;
		}
		#pragma omp critical
		{	
			outFile << kmer << '\t' << meanStd.first << '\t' << meanStd.second << '\t' << fitParameters[0] << '\t' << fitParameters[1] << '\t' << fitParameters[2] << '\t' << fitParameters[3] << '\t' << fitParameters[4] << '\t' << fitParameters[5] << '\t' << (importedEvents[i]).size() << "\t" << filteredEvents.size() << std::endl;
			pb_fit.displayProgress( prog, failed, 0 );
		}
		prog++;
	}
	outFile.close();
	std::cout << std::endl << "Done." << std::endl;

	return 0;
}
