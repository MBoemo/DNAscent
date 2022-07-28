//----------------------------------------------------------
// Copyright 2019-2020 University of Oxford
// Written by Michael A. Boemo (mb915@cam.ac.uk)
// This software is licensed under GPL-3.0.  You should have
// received a copy of the license with this software.  If
// not, please Email the author.
//----------------------------------------------------------

#ifndef DATA_IO_H
#define DATA_IO_H


#include <string>
#include <map>
#include <vector>
#include <iostream>
#include <cassert>

struct PoreParameters {

	double shift;
	double drift = 0.0;
	double scale;
	double var = 1.0;
	double eventsPerBase = 0.0;
};

class EventAlignment{

	double avg_log_emission = 0.0;
	bool spanned = false, set = false;
	unsigned int maxGap = 0;
	public:
		void recordQCs(double avg_log_emission, bool spanned, unsigned int maxGap){

			this -> avg_log_emission = avg_log_emission;
			this -> spanned = spanned;
			this -> maxGap = maxGap;
			set = true;
		}
		void printQCs(void){
			assert(set);
			std::cerr << "avg_log_emission" << " " << avg_log_emission << std::endl;
			std::cerr << "spanned" << " " << spanned << std::endl;
			std::cerr << "maxGap" << " " << maxGap << std::endl;
		}
};

struct read{

	std::string basecall, referenceSeqMappedTo, referenceMappedTo, filename, readID;
	PoreParameters scalings;
	std::vector< double > raw, normalisedEvents, eventLengths;
	std::map< unsigned int, unsigned int > refToQuery;
	std::vector< std::pair< unsigned int, unsigned int > > eventAlignment;
	std::map< unsigned int, std::pair< unsigned int, unsigned int >> eventIdx2rawIdx;
	std::map<unsigned int, double> posToScore;
	EventAlignment alignmentQCs;
	int refStart, refEnd;
	bool isReverse = false;
	public:
		void printScalings(void){

			std::cerr << "shift" << " " << scalings.shift << std::endl;
			std::cerr << "drift" << " " << scalings.drift << std::endl;
			std::cerr << "scale" << " " << scalings.scale << std::endl;
			std::cerr << "var" << " " << scalings.var << std::endl;
		}
};


/*function prototypes */
std::map< std::string, std::string > import_reference( std::string );
std::map< std::string, std::string > import_reference_pfasta( std::string );
std::vector< std::pair< double, double > > import_poreModel( std::string );
std::string getExePath(void);
std::string getGitCommit(void);
std::string writeDetectHeader(std::string, std::string, std::string, int, bool, unsigned int, unsigned int, double, bool);
std::string writeRegionsHeader(std::string, double, bool, unsigned int, unsigned int, double, double);
unsigned int sixMer2index(std::string &);


#endif
