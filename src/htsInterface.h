//----------------------------------------------------------
// Copyright 2020 University of Cambridge
// This software is licensed under GPL-3.0.  You should have
// received a copy of the license with this software.  If
// not, please Email the author.
//----------------------------------------------------------

#ifndef SRC_HTSINTERFACE_H_
#define SRC_HTSINTERFACE_H_

#include <string>
#include <utility>
#include <vector>
#include <map>
#include "../htslib/htslib/hts.h"
#include "../htslib/htslib/sam.h"

void countRecords( htsFile *, hts_idx_t *, bam_hdr_t *, int &, int , int  );
void parseCigar(bam1_t *, std::map< unsigned int, unsigned int > &, std::map< unsigned int, unsigned int > &, int &, int & );
std::string getQuerySequence( bam1_t * );
void getRefEnd(bam1_t *, int &, int & );
bool indelFastFail(bam1_t *, int, int );
std::vector<int>  ref2indels(bam1_t *, int &, int & );

#endif /* SRC_HTSINTERFACE_H_ */
