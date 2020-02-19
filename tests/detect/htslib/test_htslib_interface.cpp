//----------------------------------------------------------
// Copyright 2020 University of Cambridge
// Written by Michael A. Boemo (mb915@cam.ac.uk)
// This software is licensed under GPL-2.0.  You should have
// received a copy of the license with this software.  If
// not, please Email the author.
//----------------------------------------------------------

//NOTES: Tests DNAscent detect functions that interface with htslib

#include "../../../htslib/htslib/hts.h"
#include "../../../htslib/htslib/sam.h"
#include "../../../src/data_IO.h"
#include "../../../src/error_handling.h"
#include "../../../src/common.h" //reverse complement
#include <map>
#include <iostream>
#include <exception>

std::string getQuerySequence( bam1_t *record ){

	std::string seq;
	uint8_t *a_seq = bam_get_seq(record);
	for ( int i = 0; i < record -> core.l_qseq; i++){
		int seqInBase = bam_seqi(a_seq,i);

		switch (seqInBase) {

			case 1: seq += "A"; break;
			case 2: seq += "C"; break;
			case 4: seq += "G"; break;
			case 8: seq += "T"; break;
			case 15: seq += "N"; break;
			default: throw ParsingError();
		}
	}
	return seq;
}

void parseCigar(bam1_t *record, std::map< unsigned int, unsigned int > &ref2query, int &refStart, int &refEnd ){

	//initialise reference and query coordinates for the first match
	refStart = record -> core.pos;
	int queryPosition = 0;
	int refPosition = 0;

	const uint32_t *cigar = bam_get_cigar(record);

	if ( bam_is_rev(record) ){

		for ( int i = record -> core.n_cigar - 1; i >= 0; i--){

			const int op = bam_cigar_op(cigar[i]); //cigar operation
			const int ol = bam_cigar_oplen(cigar[i]); //number of consecutive operations

			//for a match, advance both reference and query together
			if (op == BAM_CMATCH or op == BAM_CEQUAL or op == BAM_CDIFF){

				for ( int j = refPosition; j < refPosition + ol; j++ ){

					ref2query[j] = queryPosition;
					queryPosition++;
				}
				refPosition += ol;
			}
			//for a deletion, advance only the reference position
			else if (op == BAM_CDEL or op == BAM_CREF_SKIP){

				for ( int j = refPosition; j < refPosition + ol; j++ ){

					ref2query[j] = queryPosition;
				}
				refPosition += ol;
			}
			//for insertions or soft clipping, advance only the query position
			else if (op == BAM_CSOFT_CLIP or op == BAM_CINS){

				for ( int j = refPosition; j < refPosition + ol; j++ ){

					ref2query[j] = queryPosition;
					queryPosition++;
				}
			}
			//N.B. hard clipping advances neither reference nor query, so ignore it
		}
	}
	else {

		for ( unsigned int i = 0; i < record -> core.n_cigar; ++i){

			const int op = bam_cigar_op(cigar[i]); //cigar operation
			const int ol = bam_cigar_oplen(cigar[i]); //number of consecutive operations

			//for a match, advance both reference and query together
			if (op == BAM_CMATCH or op == BAM_CEQUAL or op == BAM_CDIFF){

				for ( int j = refPosition; j < refPosition + ol; j++ ){

					ref2query[j] = queryPosition;
					queryPosition++;
				}
				refPosition += ol;
			}
			//for a deletion, advance only the reference position
			else if (op == BAM_CDEL or op == BAM_CREF_SKIP){

				for ( int j = refPosition; j < refPosition + ol; j++ ){

					ref2query[j] = queryPosition;
				}
				refPosition += ol;
			}
			//for insertions or soft clipping, advance only the query position
			else if (op == BAM_CSOFT_CLIP or op == BAM_CINS){

				for ( int j = refPosition; j < refPosition + ol; j++ ){

					ref2query[j] = queryPosition;
					queryPosition++;
				}
			}
			//N.B. hard clipping advances neither reference nor query, so ignore it
		}
	}
	refEnd = refStart + refPosition;
}

void getRefEnd(bam1_t *record, int &refStart, int &refEnd ){

	//initialise reference coordinates for the first match
	refStart = record -> core.pos;
	int refPosition = 0;

	const uint32_t *cigar = bam_get_cigar(record);

	if ( bam_is_rev(record) ){

		for ( int i = record -> core.n_cigar - 1; i >= 0; i--){

			const int op = bam_cigar_op(cigar[i]); //cigar operation
			const int ol = bam_cigar_oplen(cigar[i]); //number of consecutive operations

			//for a match
			if (op == BAM_CMATCH or op == BAM_CEQUAL or op == BAM_CDIFF){

				refPosition += ol;
			}
			//for a deletion
			else if (op == BAM_CDEL or op == BAM_CREF_SKIP){

				refPosition += ol;
			}
			//for insertions, advance only the query position so skip
			//N.B. hard clipping advances neither reference nor query, so ignore it
		}
	}
	else {

		for ( unsigned int i = 0; i < record -> core.n_cigar; ++i){

			const int op = bam_cigar_op(cigar[i]); //cigar operation
			const int ol = bam_cigar_oplen(cigar[i]); //number of consecutive operations

			//for a match, advance both reference and query together
			if (op == BAM_CMATCH or op == BAM_CEQUAL or op == BAM_CDIFF){

				refPosition += ol;
			}
			//for a deletion, advance only the reference position
			else if (op == BAM_CDEL or op == BAM_CREF_SKIP){

				refPosition += ol;
			}
			//for insertions, advance only the query position so skip
			//N.B. hard clipping advances neither reference nor query, so ignore it
		}
	}
	refEnd = refStart + refPosition;
}

int main(void){

	std::string fn_bam = "alignments.sorted.bam";
	std::string fn_reference = "reference.fasta";

	htsFile* bam_fh;
	hts_idx_t* bam_idx;
	bam_hdr_t* bam_hdr;
	hts_itr_t* itr;

	//import fasta reference
	std::map< std::string, std::string > reference = import_reference_pfasta( fn_reference );

	//load the bam
	std::cout << "Opening bam file... ";
	bam_fh = sam_open((fn_bam).c_str(), "r");
	if (bam_fh == NULL) throw IOerror(fn_bam);

	//load the index
	bam_idx = sam_index_load(bam_fh, (fn_bam).c_str());
	if (bam_idx == NULL) throw IOerror("index for "+fn_bam);

	//load the header
	bam_hdr = sam_hdr_read(bam_fh);
	std::cout << "ok." << std::endl;

	//build an iterator for all reads in the bam file
	const char *allReads = ".";
	itr = sam_itr_querys(bam_idx,bam_hdr,allReads);
	int result;

	do {
		//initialise the record and get the record from the file iterator
		bam1_t *record = bam_init1();
		result = sam_itr_next(bam_fh, itr, record);

		//get the read name (which will be the ONT readID from Albacore basecall)
		const char *queryName = bam_get_qname(record);
		if (queryName == NULL) continue;
		std::string s_queryName(queryName);
		std::cout <<"=================================================================" << std::endl;
		std::cout <<"Read name: " << s_queryName << std::endl;

		//get the name of the reference mapped to
		std::string mappedTo(bam_hdr -> target_name[record -> core.tid]);

		std::cout << "Mapped to reference: " << mappedTo << std::endl;

		int mappingQual = record -> core.qual;
		int refStart,refEnd;
		getRefEnd(record,refStart,refEnd);

		std::cout <<"Mapping quality: " << mappingQual << std::endl;
		std::cout <<"Mapped between reference coords: " << refStart << ", " << refEnd << std::endl;


		//fetch the basecall from the bam file
		std::string basecall = getQuerySequence(record);

		/*get the subsequence of the reference this read mapped to */
		std::string referenceSeqMappedTo = reference.at(mappedTo).substr(refStart, refEnd - refStart);

		//account for reverse complements
		if ( bam_is_rev(record) ){

			basecall = reverseComplement( basecall );
			referenceSeqMappedTo = reverseComplement( referenceSeqMappedTo );
		}

		//iterate on the cigar string to fill up the reference-to-query coordinate map
		std::map<unsigned int,unsigned int> refToQuery;
		parseCigar(record, refToQuery, refStart, refEnd);

		//print ref2query results
		for (auto dic = refToQuery.begin(); dic != refToQuery.end(); dic++){

			//print the position alignment
			std::cout << dic -> first << " --> " << dic -> second << std::endl;

			//print the bases
			std::string refBase, queryBase;
			if (dic -> first > referenceSeqMappedTo.length()-1) refBase = "X";
			else refBase = referenceSeqMappedTo.substr(dic -> first,1);
			if (dic -> second > basecall.length()-1) refBase = "X";
			else queryBase = basecall.substr(dic -> second,1);
			std::cout << refBase << " --> " << queryBase << std::endl;
		}

		bam_destroy1(record);

	} while(result > 0);
}
