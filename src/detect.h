//----------------------------------------------------------
// Copyright 2019 University of Oxford
// This software is licensed under GPL-3.0.  You should have
// received a copy of the license with this software.  If
// not, please Email the author.
//----------------------------------------------------------

#ifndef DETECT_H
#define DETECT_H

#include <utility>
#include <string>
#include <vector>
#include "data_IO.h"
#include "alignment.h"
#include "tensor.h"

struct HMMdetection{

	public:
		std::map<unsigned int, std::pair<double,double>> refposToLikelihood;
		std::string stdout;
};


struct DNNdetection{

	public:
		std::map<unsigned int, std::pair<double,double>> refposToProbability;
		std::string stdout;
};


/*function prototypes */
int detect_main( int argc, char** argv );

std::vector< unsigned int > getPOIs( std::string &, int );
double sequenceProbability( std::vector <double> &, std::string &, size_t, bool, PoreParameters, size_t, size_t );
DNNdetection runCNN(std::shared_ptr<AlignedRead> , std::shared_ptr<ModelSession> , std::vector<TF_Output> );
HMMdetection llAcrossRead( read &, unsigned int );

#endif
