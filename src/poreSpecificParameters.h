//----------------------------------------------------------
// Copyright 2017 University of Oxford
// Written by Michael A. Boemo (michael.boemo@path.ox.ac.uk)
// This software is licensed under GPL-3.0.  You should have
// received a copy of the license with this software.  If
// not, please Email the author.
//----------------------------------------------------------

#ifndef PORESPECIFICPARAMETERS_H
#define PORESPECIFICPARAMETERS_H

//Initial transitions within modules (internal transitions)
double internalSS2M1 = 0.97;
double internalSS2M2 = 0.03;
double internalD2I = 0.01;
double internalI2I = 0.50;
double internalI2SS = 0.49;
double internalM12M1 = 0.51;
double internalM12SE = 0.49;
double internalM22M2 = 0.97;
double internalM22SE = 0.03;
double internalSE2I = 0.01;

//Initial transitions between modules (external transitions)
double externalD2D = 0.85;
double externalD2SS = 0.14;
double externalI2SS = 0.01;
double externalSE2D = 0.12;
double externalSE2SS = 0.87;

#endif
