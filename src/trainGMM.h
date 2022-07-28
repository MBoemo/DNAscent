//----------------------------------------------------------
// Copyright 2019 University of Oxford
// Written by Michael A. Boemo (mb915@cam.ac.uk)
// This software is licensed under GPL-3.0.  You should have
// received a copy of the license with this software.  If
// not, please Email the author.
//----------------------------------------------------------

#ifndef TRAIN_H
#define TRAIN_H

/*function prototypes */
int train_main( int argc, char** argv );
std::vector< double > gaussianMixtureEM( double, double, double, double, double, std::vector< double > &, double, int );
std::map<int,int> DBSCAN( std::vector< double > &, double, unsigned int );

#endif
