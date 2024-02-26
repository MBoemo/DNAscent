//----------------------------------------------------------
// Copyright 2019 University of Oxford
// Written by Michael A. Boemo (mb915@cam.ac.uk)
// This software is licensed under GPL-3.0.  You should have
// received a copy of the license with this software.  If
// not, please Email the author.
//----------------------------------------------------------


#include <err.h>
#include <fcntl.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>
#include <libgen.h>
#include <iostream>
#include <ctime>
#include <cmath>
#include <random>
#include <algorithm>
#include <fstream>
#include <sstream>
#include "data_IO.h"
#include "pfasta/pfasta.h"
#include "gitcommit.h"
#include "common.h"
#include "softwarepath.h"
#include "forkSense.h"
#include "config.h"


std::string writeDetectHeader(std::string alignmentFilename,
		                std::string refFilename,
				std::string indexFn,
				int threads,
				bool useHMM,
				unsigned int quality,
				unsigned int length,
				bool useGPU){

	std::string detMode = "CNN";

	std::string compMode;
	if (useGPU) compMode = "GPU";
	else compMode = "CPU";

	auto t = std::time(nullptr);
	auto tm = *std::localtime(&t);
	std::ostringstream oss;
	oss << std::put_time(&tm, "%d/%m/%Y %H:%M:%S");
	auto str = oss.str();

	std::string out;
	out += "#Alignment " + alignmentFilename + "\n";
	out += "#Genome " + refFilename + "\n";
	out += "#Index " + indexFn + "\n";
	out += "#Threads " + std::to_string(threads) + "\n";
	out += "#Compute " + compMode + "\n";
	out += "#Mode " + detMode + "\n";
	out += "#MappingQuality " + std::to_string(quality) + "\n";
	out += "#MappingLength " + std::to_string(length) + "\n";
	out += "#SystemStartTime " + str + "\n";
	out += "#Software " + std::string(executablePath) + "\n";
	out += "#Version " + std::string(VERSION) + "\n";
	out += "#Commit " + std::string(gitcommit) + "\n";

	return out;
}

std::string writeRegionsHeader(std::string detectFile,
		                double threshold,
						bool useHMM,
						unsigned int cooldown,
						unsigned int resolution,
						double probability,
						double zscore){

	std::string detMode;
	if (useHMM) detMode = "HMM";
	else detMode = "CNN";

	auto t = std::time(nullptr);
	auto tm = *std::localtime(&t);
	std::ostringstream oss;
	oss << std::put_time(&tm, "%d/%m/%Y %H:%M:%S");
	auto str = oss.str();

	std::string out;
	out += "#DetectFile " + detectFile + "\n";
	out += "#Mode " + detMode + "\n";
	out += "#CallThreshold " + std::to_string(threshold) + "\n";
	out += "#Cooldown " + std::to_string(cooldown) + "\n";
	out += "#Resolution " + std::to_string(resolution) + "\n";
	out += "#Probability " + std::to_string(probability) + "\n";
	out += "#ZScore " + std::to_string(zscore) + "\n";
	out += "#SystemStartTime " + str + "\n";
	out += "#Software " + std::string(executablePath) + "\n";
	out += "#Version " + std::string(VERSION) + "\n";
	out += "#Commit " + std::string(gitcommit) + "\n";

	return out;
}


std::map< std::string, std::string > import_reference( std::string fastaFilePath ){
	
	std::ifstream file( fastaFilePath );
	if ( not file.is_open() ) throw IOerror( fastaFilePath );

	std::string line, currentRefName;
	std::map< std::string, std::string > reference;
	int referencesInFile = 0;
	
	/*while we have a line to read in the reference file... */
	while ( std::getline( file, line ) ){

		/*if this line is a fasta header line */
		if ( line[0] == '>' ){

			currentRefName = line.substr(1);
			if ( currentRefName.find(' ') != std::string::npos ) currentRefName = currentRefName.substr(0, currentRefName.find(' '));
			reference[currentRefName] = "";
			referencesInFile++;
		}
		else {

			std::transform( line.begin(), line.end(), line.begin(), toupper );
			/*grammar check: reference should only have A,T,G,C,N */
			for ( auto it = line.begin(); it < line.end(); it++ ){

				/*ignire carriage returns */
				if ( *it == '\r' ) continue;

				if ( *it != 'A' and *it != 'T' and *it != 'G' and *it != 'C' and *it != 'N' and *it != 'U' and *it != 'R' and *it != 'Y' and *it != 'K' and *it != 'M' and *it != 'S' and *it != 'W' and *it != 'B' and *it != 'D' and *it != 'H' and *it != 'V' ){
					std::cout << "Exiting with error.  Illegal character in reference file: " << *it << std::endl;
					exit( EXIT_FAILURE );
				}
			}
			reference[currentRefName] += line;
		}
	}

	/*some warning messages for silly inputs */
	if ( referencesInFile == 0 ){
		std::cout << "Exiting with error.  No fasta header (>) found in reference file." << std::endl;
		exit( EXIT_FAILURE );
	}
	return reference;
}


std::map< std::string, std::string > import_reference_pfasta( std::string fastaFilePath ){

	std::cout << "Importing reference... ";
	std::map< std::string, std::string > reference;

	int file_descriptor =
	    strcmp(fastaFilePath.c_str(), "-") == 0 ? STDIN_FILENO : open(fastaFilePath.c_str(), O_RDONLY);
	if (file_descriptor < 0) err(1, "%s", fastaFilePath.c_str());

	pfasta_file pf;
	int l = pfasta_parse(&pf, file_descriptor);
	if (l != 0) {
		errx(1, "%s: %s", fastaFilePath.c_str(), pfasta_strerror(&pf));
	}

	pfasta_seq ps;
	while ((l = pfasta_read(&pf, &ps)) == 0) {

		std::string chromosomeName(ps.name);
		std::string chromosomeSeq(ps.seq);
		std::transform( chromosomeSeq.begin(), chromosomeSeq.end(), chromosomeSeq.begin(), toupper );
		reference[ chromosomeName ] = chromosomeSeq;
		pfasta_seq_free(&ps);
	}

	if (l < 0) {
		errx(1, "%s: %s", fastaFilePath.c_str(), pfasta_strerror(&pf));
	}

	pfasta_free(&pf);
	close(file_descriptor);

	std::cout << "ok." << std::endl;
	return reference;
}


std::string getExePath(void){

	std::string s(executablePath);
	return s;
}


std::string getGitCommit(void){

	std::string s(gitcommit);
	return s;
}


unsigned int kmer2index(std::string &kmer, unsigned int kmer_len){

	std::map<std::string, int> base2index = {{"A",0}, {"T",1}, {"G",2}, {"C",3}};

	unsigned int p = 1;
	unsigned int r = 0;
	for (size_t i = 0; i < kmer_len; i++){

		r += base2index[kmer.substr(kmer_len-i-1,1)] * p;
		p *= 4;
	}
	return r;
}


std::vector< std::pair< double, double > > import_poreModel_staticStdv( std::string poreModelFilename, unsigned int kmer_len ){

	std::string pathExe = getExePath();
	std::string modelPath = pathExe + "/pore_models/" + poreModelFilename;

	/*map that sends a 5mer or 6mer to the characteristic mean and standard deviation (a pair) */
	std::map< std::string, std::pair< double, double > > kmer2MeanStd;

	/*file handle, and delimiter between columns (a \t character in the case of ONT model files) */
	std::ifstream file( modelPath );
	if ( not file.is_open() ) throw IOerror( modelPath );

	std::string line, key, mean, std;
	std::string delim = "\t";

	/*while there's a line to read */
	while ( std::getline( file, line ) ){

		/*and the line isn't part of the header */
		if ( line[0] != '#' ){ 

			/*the kmer, mean, and standard deviation are the first, second, and third columns, respectively. */
			/*take the line up to the delimiter (\t), erase that bit, and then move to the next one */
			key = line.substr( 0, line.find( delim ) );
			line.erase( 0, line.find( delim ) + delim.length() );

			mean = line.substr( 0, line.find( "\n" ) );

			/*key the map by the kmer, and convert the mean and std strings to doubles */
			kmer2MeanStd[ key ] = std::make_pair( atof(mean.c_str()), 0.14 );						
		}
	}

	std::vector< std::pair< double, double > > indexedPoreModel(pow(4,kmer_len), std::make_pair(0,0));
	
	for(auto it = kmer2MeanStd.cbegin(); it != kmer2MeanStd.cend(); ++it){

		std::string kmer = it -> first;
		kmer2MeanStd[kmer] = std::make_pair(it->second.first, it->second.second);

		indexedPoreModel[kmer2index(kmer, kmer_len)] = std::make_pair(it->second.first, it->second.second);

	    //std::cout << it->first << " " << it->second.first << " " << it->second.second << "\n";
	}

	return indexedPoreModel;
}


std::vector< std::pair< double, double > > import_poreModel_fitStdv( std::string poreModelFilename, unsigned int kmer_len ){

	std::string pathExe = getExePath();
	std::string modelPath = pathExe + "/pore_models/" + poreModelFilename;

	/*map that sends a 5mer or 6mer to the characteristic mean and standard deviation (a pair) */
	std::map< std::string, std::pair< double, double > > kmer2MeanStd;

	/*file handle, and delimiter between columns (a \t character in the case of ONT model files) */
	std::ifstream file( modelPath );
	if ( not file.is_open() ) throw IOerror( modelPath );

	std::string line, key, mean, std;
	std::string delim = "\t";

	/*while there's a line to read */
	while ( std::getline( file, line ) ){

		/*and the line isn't part of the header */
		if ( line[0] != '#' ){ 

			/*the kmer, mean, and standard deviation are the first, second, and third columns, respectively. */
			/*take the line up to the delimiter (\t), erase that bit, and then move to the next one */
			key = line.substr( 0, line.find( delim ) );
			line.erase( 0, line.find( delim ) + delim.length() );

			mean = line.substr( 0, line.find( delim ) );
			line.erase( 0, line.find( delim ) + delim.length() );

			std = line.substr( 0, line.find( "\n" ) );

			/*key the map by the kmer, and convert the mean and std strings to doubles */
			kmer2MeanStd[ key ] = std::make_pair( atof( mean.c_str() ), atof( std.c_str() ) );
		}
	}

	std::vector< std::pair< double, double > > indexedPoreModel(pow(4,kmer_len), std::make_pair(0,0));
	
	for(auto it = kmer2MeanStd.cbegin(); it != kmer2MeanStd.cend(); ++it){

		std::string kmer = it -> first;
		kmer2MeanStd[kmer] = std::make_pair(it->second.first, it->second.second);

		indexedPoreModel[kmer2index(kmer, kmer_len)] = std::make_pair(it->second.first, it->second.second);

	    //std::cout << it->first << " " << it->second.first << " " << it->second.second << "\n";
	}

	return indexedPoreModel;
}

void parseIndex( std::string indexFilename, std::map< std::string, std::string > &readID2path ){

	std::cout << "Loading DNAscent index... ";
	std::ifstream indexFile( indexFilename );
	if ( not indexFile.is_open() ) throw IOerror( indexFilename );
	std::string line;

	//get the readID to path map
	while ( std::getline( indexFile, line) ){

		std::string readID = line.substr(0, line.find('\t'));
		std::string path = line.substr(line.find('\t')+1);
		readID2path[readID] = path;
	}
	std::cout << "ok." << std::endl;
}
