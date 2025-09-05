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
#include <type_traits>


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


std::string trim(const std::string& str) {
	size_t first = str.find_first_not_of(" \t\n\r");
    if (first == std::string::npos) return "";
    size_t last = str.find_last_not_of(" \t\n\r");
    return str.substr(first, last - first + 1);
}

std::vector<std::string> split(const std::string& s, char delim ) {

    std::vector<std::string> splitString;
    std::string token;
    std::istringstream iss(s);

    if (std::isspace(delim)) {
        while (iss >> token) {
            splitString.push_back(token);
        }
    } else {
        while (std::getline(iss, token, delim)) {
            splitString.push_back(trim(token));  // optional trim
        }
    }
    return splitString;
}




const char *get_ext(const char *filename){

	const char *ext = strrchr(filename, '.');
	if(!ext || ext == filename) return "";
	return ext + 1;
}


std::string strip_extension(const std::string& filename) {
    size_t dot_position = filename.find_last_of('.');
    if (dot_position == std::string::npos) {
        return filename;
    } else {
        return filename.substr(0, dot_position);
    }
}