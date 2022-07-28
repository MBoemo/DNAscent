//----------------------------------------------------------
// Copyright 2019 University of Oxford
// Written by Michael A. Boemo (mb915@cam.ac.uk)
// This software is licensed under GPL-3.0.  You should have
// received a copy of the license with this software.  If
// not, please Email the author.
//----------------------------------------------------------


#include <random>
#include <limits>
#include "math.h"
#include "common.h"
#include "error_handling.h"
#include <math.h>


int show_version( int, char** ){

	std::cout << "Version: " << VERSION << std::endl;
	std::cout << "Written by Michael Boemo, Department of Pathology, University of Cambridge." << std::endl;
	std::cout << "Please submit bug reports to GitHub Issues." << std::endl;
	return 0;
}


void displayProgress( int current, int total ){
/*poached from https://stackoverflow.com/questions/14539867/how-to-display-a-progress-indicator-in-pure-c-c-cout-printf */

	double progress = (double) current / (double) total;
	int barWidth = 70;

	if ( progress <= 1.0 ){

		std::cout << "[";
		int pos = barWidth * progress;
		for (int i = 0; i < barWidth; ++i) {
			if (i < pos) std::cout << "=";
			else if (i == pos) std::cout << ">";
			else std::cout << " ";
		}
		std::cout << "] " << int(progress * 100.0) << " %\r";
		std::cout.flush();
	}
}


std::vector< std::string > split( std::string s, char delim ){

	std::stringstream ssString( s );
	std::vector< std::string > splitString;
	std::string entry;

	while ( std::getline( ssString, entry, delim ) ){

		splitString.push_back( entry );
	}
	return splitString;
}


double vectorMean( std::vector< double > &obs ){

	double total = 0.0;
	for ( size_t i = 0; i < obs.size(); i++ ) total += obs[i];
	return total / (double) obs.size();
}


double vectorStdv( std::vector< double > &obs, double &mean ){

	double total = 0.0;
	for ( size_t i = 0; i < obs.size(); i++ ){
		total += pow(obs[i] - mean, 2.0);
	}
	return sqrt(total / (double) obs.size());
}


double vectorSum( std::vector< double > &obs ){

	double total = 0.0;
	for ( size_t i = 0; i < obs.size(); i++ ) total += obs[i];
	return total;
}


int argMin( std::vector< double > vec ){

	double smallest = vec[0];
	int index = 0;

	for ( unsigned int i = 1; i < vec.size(); i++ ){

		if ( vec[i] < smallest ){

			smallest = vec[i];
			index = i;
		}
	}
	return index;
}


double logistic(double input, double slope, double centre){

	return 1/(1 + exp (-slope*(abs(input)-centre))); 
}


std::vector<double> movingAvgFilter( std::vector<double> &input, unsigned int filterSize){

	std::vector<double> out;
	
	for (size_t i = filterSize/2; i < input.size() - filterSize/2; i++){
	
		double sum = 0.;
		for (size_t j = i - filterSize/2; j < i + filterSize/2; j++){
			sum += input[j];
		}
		out.push_back( sum / (double) filterSize );
	}
	return out;
}


std::vector<double> movingAvgFilterLogistic( std::vector<double> &input, unsigned int filterSize){

	std::vector<double> out;
	
	for (size_t i = filterSize/2; i < input.size() - filterSize/2; i++){
	
		double sum = 0.;
		for (size_t j = i - filterSize/2; j < i + filterSize/2; j++){
			sum += input[j];
		}
		out.push_back( logistic( sum / (double) filterSize, 20, 0.2));
	}
	return out;
}


std::vector<double> normVectorSum(std::vector<double> input){

	double sum = vectorSum(input);

	std::vector<double> output;
	for (size_t i = 0; i < input.size(); i++){
	
		output.push_back(input[i]/sum);
	}
	return output;
}


int argMax(std::vector< double > vec){

	double highest = vec[0];
	int index = 0;

	for ( unsigned int i = 1; i < vec.size(); i++ ){

		if ( vec[i] > highest ){

			highest = vec[i];
			index = i;
		}
	}
	return index;
}
