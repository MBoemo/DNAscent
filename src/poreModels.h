//----------------------------------------------------------
// Copyright 2019 University of Oxford
// Written by Michael A. Boemo (michael.boemo@path.ox.ac.uk)
// This software is licensed under GPL-2.0.  You should have
// received a copy of the license with this software.  If
// not, please Email the author.
//----------------------------------------------------------


#ifndef POREMODELS_INL
#define POREMODELS_INL

#include <utility>
#include <string>
#include <map>

std::map< std::string, std::pair< double, double > > loadModel_SixMerONT(void);
std::map< std::string, std::pair< double, double > > loadModel_BrdUFull(void);
double KLdivergence( double mu1, double sigma1, double mu2, double sigma2 );
std::map< std::string, std::pair< double, double > > buildAnalogueModel(double KLthreshold, bool excludeCpG);

#endif
