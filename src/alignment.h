//----------------------------------------------------------
// Copyright 2019-2020 University of Oxford
// This software is licensed under GPL-3.0.  You should have
// received a copy of the license with this software.  If
// not, please Email the author.
//----------------------------------------------------------

#ifndef ALIGN_H
#define ALIGN_H

#include <fstream>
#include <math.h>
#include <stdlib.h>
#include <limits>
#include "../fast5/include/fast5.hpp"
#include <memory>
#include <utility>
#include "reads.h"


int align_main( int argc, char** argv );
void eventalign( DNAscent::read &, unsigned int);

#endif
