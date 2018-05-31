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
double internalM12I = 0.001;
double internalI2I = 0.001;
double internalM12M1 = 0.4;

//Initial transitions between modules (external transitions)
double externalD2D = 0.3;
double externalD2M1 = 0.7;
double externalI2M1 = 0.999;
double externalM12D = 0.0025;
double externalM12M1 = 0.5965;

#endif
