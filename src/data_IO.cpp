//----------------------------------------------------------
// Copyright 2017 University of Oxford
// Written by Michael A. Boemo (michael.boemo@path.ox.ac.uk)
// This software is licensed under GPL-3.0.  You should have
// received a copy of the license with this software.  If
// not, please Email the author.
//----------------------------------------------------------


#include <iostream>
#include <random>
#include <algorithm>
#include <fstream>
#include <sstream>
#include "data_IO.h"
#include "error_handling.h"
#include "poreModels.h"


std::string import_reference( std::string fastaFilePath ){
	
	std::ifstream file( fastaFilePath );

	if ( not file.is_open() ) throw IOerror( fastaFilePath );

	std::string line;
	std::string reference;
	int referencesInFile = 0;
	
	/*while we have a line to read in the reference file... */
	while ( std::getline( file, line ) ){

		/*if this line is a fasta header line */
		if ( line[0] == '>' ) referencesInFile++;

		/*if this line is not a fasta header line AND is the first entry */
		if ( line[0] != '>' and referencesInFile == 1 ){

			/*append the line to the reference string */
			reference.append( line );
		}
	}

	/*some warning messages for silly inputs */
	if ( referencesInFile > 1 ){
		std::cout << "WARNING: reference fasta file contains multiple entries.  Only using the first one." << std::endl;
	}
	if ( referencesInFile == 0 ){
		std::cout << "Exiting with error.  No fasta header (>) found in reference file." << std::endl;
		exit( EXIT_FAILURE );
	}

	/*put all characters in the reference to upper case */
	std::transform( reference.begin(), reference.end(), reference.begin(), toupper );

	/*grammar check: reference should only have A,T,G,C,N */
	for ( auto it = reference.begin(); it < reference.end(); it++ ){

		/*ignire carriage returns */
		if ( *it == '\r' ) continue;

		if ( *it != 'A' and *it != 'T' and *it != 'G' and *it != 'C' and *it != 'N' ){
			std::cout << "Exiting with error.  Illegal character in reference file: " << *it << std::endl;
			exit( EXIT_FAILURE );
		}
	}
	return reference;
}


std::map< std::string, std::pair< double, double > > import_poreModel( std::string poreModelFilePath ){

	/*map that sends a 5mer or 6mer to the characteristic mean and standard deviation (a pair) */
	std::map< std::string, std::pair< double, double > > kmer2MeanStd;

	/*file handle, and delimiter between columns (a \t character in the case of ONT model files) */
	std::ifstream file( poreModelFilePath );

	if ( not file.is_open() ) throw IOerror( poreModelFilePath );

	std::string line, key, mean, std;
	std::string delim = "\t";

	/*while there's a line to read */
	while ( std::getline( file, line ) ){

		/*and the line isn't part of the header */
		if ( line.substr(0,4) != "kmer" && line[0] != '#' ){ 

			/*the kmer, mean, and standard deviation are the first, second, and third columns, respectively. */
			/*take the line up to the delimiter (\t), erase that bit, and then move to the next one */
			key = line.substr( 0, line.find( delim ) );
			line.erase( 0, line.find( delim ) + delim.length() );

			mean = line.substr( 0, line.find( delim ) );
			line.erase( 0, line.find( delim ) + delim.length() );

			std = line.substr( 0, line.find( delim ) );

			/*key the map by the kmer, and convert the mean and std strings to doubles */
			kmer2MeanStd[ key ] = std::make_pair( atof( mean.c_str() ), atof( std.c_str() ) );
		}
	}

	/*if you need to print out the map for debugging purposes 
	for(auto it = kmer2MeanStd.cbegin(); it != kmer2MeanStd.cend(); ++it){
	    std::cout << it->first << " " << it->second.first << " " << it->second.second << "\n";
	}*/
	
	return kmer2MeanStd;
}


void export_poreModel( std::map< std::string, std::vector< std::pair< double, double > > > &trainedMap, std::string &outputFilename ){

	std::ofstream outFile;
	outFile.open( outputFilename );
	if ( not outFile.is_open() ) throw IOerror( outputFilename );

	/*write headers */
	outFile << "#Analogue pore model file trained by Osiris." << std::endl;
	outFile << "kmer\ttrMu\ttrSig\toriMu\toriSig" << std::endl;

	/*uniform random number generator */
	std::random_device rd{};
	std::mt19937 gen{rd()};

	/*for each 5mer that we trained in the model */
	for( auto iter = trainedMap.cbegin(); iter != trainedMap.cend(); ++iter ){

		std::vector< double > samples;
		samples.reserve( 1000*(iter -> second).size() );

		for ( auto muStdPair = (iter -> second).begin(); muStdPair< (iter -> second).end(); muStdPair++ ){

			/*generator for the normal distribution corresponding to this training value */
			std::normal_distribution<> nd{ muStdPair -> first, muStdPair -> second };
			
			/*draw 1000 samples */
			for ( int i = 0; i < 1000; i++ ){

				samples.push_back( nd(gen) );
			}
		}

		/*use the samples to calculate a consensus mean and stdv between trained values */
		double sum = 0.0;
		for ( int i = 0; i < samples.size(); i++ ){

			sum += samples[i];
		}
		double mean = sum / (double) samples.size();

		double stdvSum = 0.0;
		for ( int i = 0; i < samples.size(); i++ ){

			stdvSum += pow(samples[i] - mean,2);
		}
		double stdv = sqrt( stdvSum / (double) samples.size() );
		std::string fiveMer = (iter -> first);
		fiveMer.replace( (iter -> first).find('B'), 1, "T" );

		outFile << iter -> first << "\t" << mean << "\t" << stdv << "\t" << FiveMer_model[fiveMer].first << "\t" << FiveMer_model[fiveMer].second << std::endl;
	}
	outFile.close();
}


std::pair< std::string, std::vector< read > > getTrainingFrom_foh( std::string &trainingGroup ){
	
	std::string line, sevenMerName, event, bound;
	std::stringstream trainingGroupStream;

	read currentRead;
	std::vector< read > trainingGroupReads;

	/*read this training group from the foh file into a stream */
	trainingGroupStream.str( trainingGroup );

	/*get the 7mer name */
	std::getline( trainingGroupStream, line );
	if ( line == "" ) std::getline( trainingGroupStream, line );
	assert( line.substr(0,1) == ">" );
	sevenMerName = line.erase( 0, 1 );
	
	/*go through each read */
	while ( std::getline( trainingGroupStream, line ) ){

		/*the basecall line */
		currentRead.basecalls = line;

		/*the reference bounds line */
		std::getline( trainingGroupStream, line );
		(currentRead.bounds_reference).first = atoi( (line.substr( 0, line.find(' ') )).c_str() );
		(currentRead.bounds_reference).second = atoi( (line.substr( line.find(' ') + 1, line.size() - line.find(' ') )).c_str() );

		/*the query bounds line */
		std::getline( trainingGroupStream, line );
		(currentRead.bounds_query).first = atoi( (line.substr( 0, line.find(' ') )).c_str() );
		(currentRead.bounds_query).second = atoi( (line.substr( line.find(' ') + 1, line.size() - line.find(' ') )).c_str() );

		/*the raw signal line */
		std::getline( trainingGroupStream, line );
		std::vector< double > rawSignals;
		std::istringstream ss( line );
		std::string event;
		while ( std::getline( ss, event, ' ' ) ){

			rawSignals.push_back( atof( event.c_str() ) );
		}
		currentRead.raw = rawSignals;
		trainingGroupReads.push_back( currentRead );
	}
	return std::make_pair( sevenMerName, trainingGroupReads );
}


read getDetectionFrom_fdh( std::string &detectionGroup ){
	
	std::stringstream detectionGroupStream;	
	read currentRead;
	std::string line;

	/*read this training group from the foh file into a stream */
	detectionGroupStream.str( detectionGroup );

	/*get the filename */
	std::getline( detectionGroupStream, line );
	if ( line == "" ) std::getline( detectionGroupStream, line );
	assert( line.substr(0,1) == ">" );
	currentRead.filename = line.erase( 0, 1 );

	/*get the basecall */
	std::getline( detectionGroupStream, line );
	currentRead.basecalls = line;

	/*extract events */
	std::getline( detectionGroupStream, line );
	std::istringstream ss( line );
	std::string event;	
	while ( std::getline( ss, event, ' ' ) ){

		(currentRead.raw).push_back( atof( event.c_str() ) );
	}

	return currentRead;
}
