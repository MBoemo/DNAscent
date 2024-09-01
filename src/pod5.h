//----------------------------------------------------------
// Copyright 2024 University of Cambridge
// This software is licensed under GPL-3.0.  You should have
// received a copy of the license with this software.  If
// not, please Email the author.
//----------------------------------------------------------

#ifndef POD5_H
#define POD5_H

#include <vector>
#include <string>
#include <unordered_map>
#include <omp.h>
#include "reads.h"

void pod5_getSignal(DNAscent::read &);
void pod5_getSignal_batch(std::vector<DNAscent::read *>);
std::vector< std::string > pod5_extract_readIDs(std::string);

#endif
