//----------------------------------------------------------
// Copyright 2019 University of Oxford
// This software is licensed under GPL-3.0.  You should have
// received a copy of the license with this software.  If
// not, please Email the author.
//----------------------------------------------------------


#ifndef COMMON_H
#define COMMON_H

#define VERSION "4.0.1"

#include <algorithm>
#include <vector>
#include <utility>
#include <iterator>
#include <map>
#include <sstream>
#include <iostream>
#include <chrono>
#include <iomanip>

int show_version( int, char** );

class progressBar{

	private:
		std::chrono::time_point<std::chrono::steady_clock> _startTime;
		std::chrono::time_point<std::chrono::steady_clock> _currentTime;
		unsigned int maxNumber;
		unsigned int barWidth = 35;
		unsigned int _digits;
		bool _withFail;

	public:
		progressBar( unsigned int maxNumber, bool withFail ){
		
			this -> maxNumber = maxNumber;
			_digits = std::to_string( maxNumber).length() + 1;
			_startTime = std::chrono::steady_clock::now();
			_withFail = withFail;
		}
		void displayProgress( unsigned int currentNumber, unsigned int failed, int failedEvents ){

			_currentTime = std::chrono::steady_clock::now();
			 std::chrono::duration<double> elapsedTime = _currentTime - _startTime;

			double progress = (double) currentNumber / (double) maxNumber;
			
			if ( progress <= 1.0 ){

				std::cout << "[";
				unsigned int pos = barWidth * progress;
				for (unsigned int i = 0; i < barWidth; ++i) {
					if (i < pos) std::cout << "=";
					else if (i == pos) std::cout << ">";
					else std::cout << " ";
				}
				std::cout << "] " << std::right << std::setw(3) << int(progress * 100.0) << "%  ";

				std::cout << std::right << std::setw(_digits) << currentNumber << "/" << maxNumber << "  ";


				unsigned int estTimeLeft = elapsedTime.count() * ( (double) maxNumber / (double) currentNumber - 1.0 );			
				unsigned int hours = estTimeLeft / 3600;
				unsigned int mins = (estTimeLeft % 3600) / 60;
				unsigned int secs = (estTimeLeft % 3600) % 60;

				if (_withFail){

					std::cout << std::right << std::setw(2) << hours << "hr" << std::setw(2) << mins << "min" << std::setw(2) << secs << "sec  ";
					std::cout << "failed: " << std::right << std::setw(_digits) << failed << std::setw(3) << "\r";

					//testing
					//std::cout << "  fe: " << std::right << std::setw(_digits) << failedEvents << std::setw(3) << "\r";
				}
				else{

					std::cout << std::right << std::setw(2) << hours << "hr" << std::setw(2) << mins << "min" << std::setw(2) << secs << "sec  " << "\r";
				}
				std::cout.flush();
			}
		} 
};


inline std::string reverseComplement( std::string DNAseq ){

	std::reverse( DNAseq.begin(), DNAseq.end() );
	std::string revComp;

	for ( std::string::iterator i = DNAseq.begin(); i < DNAseq.end(); i++ ){
	
		switch( *i ){
			case 'A' :
				revComp += 'T';
				break;
			case 'T' :
				revComp += 'A';
				break;
			case 'G' :
				revComp += 'C';
				break;
			case 'C' :
				revComp += 'G';
				break;
			case 'U' :
				revComp += 'A';
				break;
			case 'Y' :
				revComp += 'R';
				break;
			case 'R' :
				revComp += 'Y';
				break;
			case 'K' :
				revComp += 'M';
				break;
			case 'M' :
				revComp += 'K';
				break;
			case 'B' :
				revComp += 'V';
				break;
			case 'D' :
				revComp += 'H';
				break;
			case 'H' :
				revComp += 'D';
				break;
			case 'V' :
				revComp += 'B';
				break;
			case 'N' :
				revComp += 'N';
				break;
			case 'W' :
				revComp += 'W';
				break;
			case 'S' :
				revComp += 'S';
				break;
			default:
				std::cout << "Exiting with error.  Invalid character " << *i << " passed to reverse complement function - must be IUPAC character." << std::endl;
				exit( EXIT_FAILURE );
		}
	}
	return revComp;
}


inline std::string complement( std::string DNAseq ){

	std::string comp;

	for ( std::string::iterator i = DNAseq.begin(); i < DNAseq.end(); i++ ){
	
		switch( *i ){
			case 'A' :
				comp += 'T';
				break;
			case 'T' :
				comp += 'A';
				break;
			case 'G' :
				comp += 'C';
				break;
			case 'C' :
				comp += 'G';
				break;
			default:
				std::cout << "Exiting with error.  Invalid character passed to reverse complement function.  Must be A, T, G, or C." << std::endl;
				exit( EXIT_FAILURE );
		}
	}
	return comp;
}


/*function prototypes */
void displayProgress( int, int );
std::vector< std::string > split( std::string, char );
int argMin( std::vector< double > );
int argMax( std::vector< double > );
double vectorMean( std::vector< double > & );
double vectorStdv( std::vector< double > &, double & );
double vectorSum( std::vector< double > & );
std::vector<double> movingAvgFilter(std::vector<double> &, unsigned int);
std::vector<double> movingAvgFilterLogistic(std::vector<double> &, unsigned int);
std::vector<double> normVectorSum(std::vector<double>);

#endif
