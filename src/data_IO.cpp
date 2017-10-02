//----------------------------------------------------------
// Copyright 2017 University of Oxford
// Written by Michael A. Boemo (michael.boemo@path.ox.ac.uk)
//----------------------------------------------------------


#include "data_IO.h"
#include "error_handling.h"


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


void export_poreModel( std::map< std::string, std::vector< double > > &trainedMap, std::string &outputFilename ){

	std::ofstream outFile;
	outFile.open( outputFilename );
	if ( not outFile.is_open() ) throw IOerror( outputFilename );

	/*write headers */
	outFile << "#Analogue pore model file trained by Osiris." << std::endl;
	outFile << "kmer\ttrMu\ttrSig\toriMu\toriSig" << std::endl;

	for( auto iter = trainedMap.cbegin(); iter != trainedMap.cend(); ++iter ){

		outFile << iter -> first << "\t" << ( iter -> second )[ 0 ] << "\t" << ( iter -> second )[ 1 ] << "\t" << ( iter -> second )[ 2 ] << "\t" << ( iter -> second )[ 3 ] << std::endl;
	}
	outFile.close();
}


std::map< std::string, std::vector< std::vector< double > > > import_foh( std::string fohFilePath ){
	
	std::cout << "Importing training data..." << std::endl;
	std::ifstream file( fohFilePath );
	
	if ( not file.is_open() ) throw IOerror( fohFilePath );

	std::string line;
	std::string reference;
	std::string delim = " ";
	std::string key;
	std::string event;

	std::map< std::string, std::vector< std::vector< double > > > kmer2normalisedReads;
	
	/*while we have a line to read in the reference file... */
	while ( std::getline( file, line ) ){

		if ( line[0] == '>' ){

			key = line.erase( 0, 1 );
		}
		/*if this line is not a fasta header line */
		else{

			std::vector< double > Events;

			std::istringstream ss( line );
			std::string event;
			while ( std::getline( ss, event, ' ' ) ){

				Events.push_back( atof( event.c_str() ) );
			}
			kmer2normalisedReads[ key ].push_back( Events );
		}
	}
	std::cout << "Done." << std::endl;
	return kmer2normalisedReads;

}


std::vector< detectionTuple > import_fdh( std::string &fdhFilePath ){
	
	std::cout << "Importing detection data..." << std::endl;
	std::ifstream file( fdhFilePath );
	
	if ( not file.is_open() ) throw IOerror( fdhFilePath );

	std::string line;
	std::string event;

	std::vector< detectionTuple > detectionTuples;
	
	/*while we have a line to read in the reference file... */
	while ( std::getline( file, line ) ){

		if ( line[0] == '>' ){

			detectionTuple readTuple;

			/*extract filename */
			readTuple.filename = line.erase( 0, 1 );

			/*extract basecalls */
			std::getline( file, line );
			readTuple.basecalls = line;

			/*extract events */
			std::getline( file, line );
			std::istringstream ss( line );
			std::string event;
			while ( std::getline( ss, event, ' ' ) ){

				(readTuple.events).push_back( atof( event.c_str() ) );
			}
			detectionTuples.push_back( readTuple );
		}
	}

	std::cout << "Done." << std::endl;
	return detectionTuples;

}
