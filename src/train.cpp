//----------------------------------------------------------
// Copyright 2019 University of Oxford
// Written by Michael A. Boemo (mb915@cam.ac.uk)
// This software is licensed under GPL-2.0.  You should have
// received a copy of the license with this software.  If
// not, please Email the author.
//----------------------------------------------------------

#include <exception>
#include <algorithm>
#include <math.h>
#include <vector>
#include <iostream>
#include <fstream>
#include "common.h"
#include "data_IO.h"
#include "error_handling.h"
#include "event_handling.h"
#include "probability.h"
#include "poreModels.h"
#include "train.h"

static const char *help=
"train: DNAscent executable that determines the mean and standard deviation of a base analogue's current.\n"
"To run DNAscent train, do:\n"
"  ./DNAscent train -d /path/to/eventalign.nanopolish -o output.model\n"
"Required arguments are:\n"
"  -d,--trainingData         path to training data from nanopolish eventalign,\n"
"  -o,--output               path to the output trained pore model file.\n"
"Optional arguments are:\n"
"  -pi,                      mixing parameter for BrdU (default is 0.5),\n"
"  -m,--max-reads            maximum number of reads to consider (default is 100000),\n"
"  -e,--max-events           maximum number of events per 6mer to consider (default is 10000),\n"
"  -t,--threads              number of threads (default is 1 thread).\n"
"Written by Michael Boemo, Department of Pathology, University of Cambridge.\n"
"Please submit bug reports to GitHub Issues (https://github.com/MBoemo/DNAscent/issues).";

struct Arguments {

	std::string trainingOutputFilename;
	bool logFile;
	std::string logFilename;
	int threads;
	std::string eventalignFilename;
	unsigned int maxReads;
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


std::map<int,int> DBSCAN( std::vector< double > events, double epsilon, unsigned int minPoints ){

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
	for ( unsigned int i = 0; i < data.size(); i++ ){

		logLikelihood_Old += eln(pi1 * normalPDF(mu1,sigma1,data[i]) + pi2 * normalPDF(mu2,sigma2,data[i])); 
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
		//pi1 = Nk1 / (double) data.size();
		//pi2 = Nk2 / (double) data.size();

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

			logLikelihood_New += eln(pi1 * normalPDF(mu1,sigma1,data[i]) + pi2 * normalPDF(mu2,sigma2,data[i])); 
		}

		improvement = logLikelihood_New - logLikelihood_Old;
		logLikelihood_Old = logLikelihood_New;
		iterations++;
		if (iterations > maxIter) break;
	}
	return { pi1, mu1, sigma1, pi2, mu2, sigma2 };
}


std::vector< double > gaussianMixtureEM( double mu1, double sigma1, double mu2, double sigma2, std::vector< double > &data, double tolerance, int maxIter ){

	std::vector< std::vector< double > > Z( 2, std::vector< double >( data.size(), 0 ) );

	double pi1 = 0.5;
	double pi2 = 0.5;
	double total1, total2, logLikelihood_Old, logLikelihood_New;

	/*INITIALISATION - calculuate the log likelihood using the paramters initially passed */
	logLikelihood_Old = 0;
	for ( unsigned int i = 0; i < data.size(); i++ ){

		logLikelihood_Old += eln(pi1 * normalPDF(mu1,sigma1,data[i]) + pi2 * normalPDF(mu2,sigma2,data[i])); 
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

			logLikelihood_New += eln(pi1 * normalPDF(mu1,sigma1,data[i]) + pi2 * normalPDF(mu2,sigma2,data[i])); 
		}

		improvement = logLikelihood_New - logLikelihood_Old;
		logLikelihood_Old = logLikelihood_New;
		iterations++;
		if (iterations > maxIter) break;
	}
	return { pi1, mu1, sigma1, pi2, mu2, sigma2 };
}


//https://stackoverflow.com/questions/843154/fastest-way-to-find-the-number-of-lines-in-a-text-c
unsigned int FileRead( std::istream &is, std::vector<char> &buff ){
	is.read( &buff[0], buff.size() );
	return is.gcount();
}

//https://stackoverflow.com/questions/843154/fastest-way-to-find-the-number-of-lines-in-a-text-c
unsigned int CountLines( const std::vector<char> &buff, int sz ){
	int newlines = 0;
	const char * p = &buff[0];
	for ( int i = 0; i < sz; i++ ){
		if ( p[i] == '\n' ){
			newlines++;
		}
	}
	return newlines;
}


int train_main( int argc, char** argv ){

	Arguments trainArgs = parseTrainingArguments( argc, argv );

	/*open output file */
	std::ofstream outFile( trainArgs.trainingOutputFilename );
	if ( not outFile.is_open() ) throw IOerror( trainArgs.trainingOutputFilename );

	/*get events from the work file */
	std::cout << "Trawling through eventalign... ";

	std::string line;
	int prog, failed;

	/*fudge for openmp */
	std::map< int, std::string > indexToSixmer;
	std::map< std::string, int > sixmerToIndex;
	int index = 0;
	for ( auto i = thymidineModel.cbegin(); i != thymidineModel.cend(); i++ ){

		indexToSixmer[index] = i -> first;
		sixmerToIndex[i->first] = index;
		index++;
	}

	std::vector< std::vector< double > > importedEvents( 4096 );

	/*
	//https://stackoverflow.com/questions/843154/fastest-way-to-find-the-number-of-lines-in-a-text-c
  	std::ifstream inFile(trainArgs.eventalignFilename); 
	const int SZ = 1024 * 1024;
        std::vector <char> buff( SZ );
        int lineCount = 0;
        while( int cc = FileRead( inFile, buff ) ) {
            lineCount += CountLines( buff, cc );
        }
	inFile.close();
	progressBar pb_read( lineCount, true );
	*/
	progressBar pb_read( trainArgs.maxReads, true );
	unsigned int readsRead = 0;
	std::string readIdx = "";

	std::ifstream eventFile(trainArgs.eventalignFilename);
	if ( not eventFile.is_open() ) throw IOerror( trainArgs.eventalignFilename );

	std::getline( eventFile, line);//throw away the header

	while ( std::getline( eventFile, line) ){

		std::istringstream ss( line );
		std::string sixMer, entry;
		double eventMean = 0.0, eventLength = 0.0;

		int col = 0;
		while ( std::getline( ss, entry, '\t' ) ){

			if ( col == 9 ){

				sixMer = entry;
				break;
			}
			else if ( col == 3 ){

				if (entry != readIdx){

					readsRead++;
					pb_read.displayProgress( readsRead, 0, 0 );
					readIdx = entry;
				}
			}
			else if ( col == 6 ){

				eventMean = atof( entry.c_str() );
			}
			else if ( col == 8 ){

				eventLength = atof( entry.c_str() );
			}
			col++;
		}

		assert (eventMean != 0.0 and eventLength != 0.0);

		if ( eventLength >= 0.002 and importedEvents[sixmerToIndex[sixMer]].size() < trainArgs.maxEvents ){

			importedEvents[sixmerToIndex[sixMer]].push_back( eventMean );
		}
		if (readsRead > trainArgs.maxReads) break;
	}
	eventFile.close();
	std::cout << "ok." << std::endl;
	std::cout << "Fitting..." << std::endl;
	
	/*fit a mixture model to the events that aligned to each position in the reference */
	outFile << "6mer" << '\t' << "ONT_mean" << '\t' << "ONT_stdv" << '\t' << "pi_1" << '\t' << "mean_1" << '\t' << "stdv_1" << '\t' << "pi_2" << '\t' << "mean_2" << '\t' << "stdv_2" << std::endl;
	progressBar pb_fit( importedEvents.size(),true );

	#pragma omp parallel for schedule(dynamic) shared(pb_fit, indexToSixmer, thymidineModel, prog, failed, outFile, importedEvents, trainArgs) num_threads(trainArgs.threads)
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

		std::string sixMer = indexToSixmer[i];
		double mu1, stdv1, mu2, stdv2;

		/*get the ONT distribution for the mixture */
		mu1 = thymidineModel[sixMer].first;
		stdv1 = thymidineModel[sixMer].second;

		/*make a second distribution that's similar to the ONT distribution */
		mu2 = thymidineModel[sixMer].first;
		stdv2 = 2*thymidineModel[sixMer].second;

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
			outFile << sixMer << '\t' << thymidineModel[sixMer].first << '\t' << thymidineModel[sixMer].second << '\t' << fitParameters[0] << '\t' << fitParameters[1] << '\t' << fitParameters[2] << '\t' << fitParameters[3] << '\t' << fitParameters[4] << '\t' << fitParameters[5] << '\t' << (importedEvents[i]).size() << "\t" << filteredEvents.size() << std::endl; 
			pb_fit.displayProgress( prog, failed, 0 );
		}
		prog++;
	}
	outFile.close();
	std::cout << std::endl << "Done." << std::endl;

	return 0;
}
