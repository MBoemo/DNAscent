//----------------------------------------------------------
// Copyright 2017 University of Oxford
// Written by Michael A. Boemo (michael.boemo@path.ox.ac.uk)
// This software is licensed under GPL-2.0.  You should have
// received a copy of the license with this software.  If
// not, please Email the author.
//----------------------------------------------------------

#ifndef ERROR_HANDLING_H
#define ERROR_HANDLING_H

#include "scrappie/event_detection.h"
#include "poreModels.h"
#include "common.h"
#include "data_IO.h"
#include <math.h>
#include <cmath>
#include <stdlib.h>
#include <assert.h>
#include <exception>
#include <string.h>

struct IOerror : public std::exception {
	std::string badFilename;	
	IOerror( std::string s ){

		badFilename = s;
	}
	const char* what () const throw () {
		const char* message = "Could not open file: ";
		const char* specifier = badFilename.c_str();
		char* result;
		result = static_cast<char*>(calloc(strlen(message)+strlen(specifier)+1, sizeof(char)));
		strcpy( result, message);
		strcat( result, specifier );

		return result;
	}	
};


struct InvalidOption : public std::exception {
	std::string badOption;	
	InvalidOption( std::string s ){

		badOption = s;
	}
	const char* what () const throw () {
		const char* message = "Invalid option passed: ";
		const char* specifier = badOption.c_str();
		char* result;
		result = static_cast<char*>(calloc(strlen(message)+strlen(specifier)+1, sizeof(char)));
		strcpy( result, message);
		strcat( result, specifier );

		return result;
	}	
};


struct InsufficientArguments : public std::exception {
	const char * what () const throw () {
		return "Insufficient number of arguments passed to executable.";
	}
};


struct BadStrandDirection : public std::exception {
	const char * what () const throw () {
		return "Unrecognised strand direction.";
	}
};


struct FastaFormatting : public std::exception {
	const char * what () const throw () {
		return "Reference file is not in the correct format - fasta format required.";
	}
};


struct BadFast5Field : public std::exception {
	const char * what () const throw () {
		return "Fast5 field could not be opened.";
	}
};


struct MismatchedDimensions : public std::exception {
	const char * what () const throw () {
		return "Gaussian elimination on A*x=b.  Rows in A must equal length of b.";
	}
};


struct ParsingError : public std::exception {
	const char * what () const throw () {
		return "Parsing error reading BAM file.";
	}
};

#endif

