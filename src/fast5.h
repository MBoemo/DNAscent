//----------------------------------------------------------
// Copyright 2019-2020 University of Oxford
// This software is licensed under GPL-3.0.  You should have
// received a copy of the license with this software.  If
// not, please Email the author.
//----------------------------------------------------------

#ifndef FAST5_H
#define FAST5_H

#include <vector>
#include <string>
#include "reads.h"

void fast5_getSignal( DNAscent::read & );
std::vector<std::string> fast5_extract_readIDs(std::string);

#endif
