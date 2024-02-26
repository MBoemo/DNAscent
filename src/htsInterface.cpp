//----------------------------------------------------------
// Copyright 2020 University of Cambridge
// Written by Michael A. Boemo (mb915@cam.ac.uk)
// This software is licensed under GPL-3.0.  You should have
// received a copy of the license with this software.  If
// not, please Email the author.
//----------------------------------------------------------

#include "htsInterface.h"
#include <iostream>
#include "error_handling.h"


void countRecords( htsFile *bam_fh, hts_idx_t *bam_idx, bam_hdr_t *bam_hdr, int &numOfRecords, int minQ, int minL ){

	std::cout << "Scanning bam file...";
	hts_itr_t* itr = sam_itr_querys(bam_idx,bam_hdr,".");
	int result;

	do {
		bam1_t *record = bam_init1();
		result = sam_itr_next(bam_fh, itr, record);
		int refStart,refEnd;
		getRefEnd(record,refStart,refEnd);
		int queryLen = record -> core.l_qseq;
		if ( (record -> core.qual >= minQ) and (refEnd - refStart >= minL) and queryLen != 0) numOfRecords++;
		bam_destroy1(record);
	} while (result > 0);

	//cleanup
	sam_itr_destroy(itr);
	std::cout << "ok." << std::endl;
}


bool indelFastFail(bam1_t *record, int maxI, int maxD ){
	//Covered in: tests/detect/htslib

	const uint32_t *cigar = bam_get_cigar(record);

	for ( unsigned int i = 0; i < record -> core.n_cigar; ++i){

		const int op = bam_cigar_op(cigar[i]); //cigar operation
		const int ol = bam_cigar_oplen(cigar[i]); //number of consecutive operations

		//deletions
		if (op == BAM_CDEL or op == BAM_CREF_SKIP){
			if (ol >= maxD){
				return true;
			}
		}
		//for insertions or soft clipping, advance only the query position
		else if (op == BAM_CINS){
			if (ol >= maxI){
				return true;
			}
		}
	}
	return false;
}


std::vector<int> ref2indels(bam1_t *record, int &refStart, int &refEnd ){
	//Covered in: tests/detect/htslib

	//initialise reference and query coordinates for the first match
	//refStart = record -> core.pos;
	int queryPosition = 0;
	int refPosition = 0;

	std::vector<int> indels(refEnd - refStart + 1, 0);

	const uint32_t *cigar = bam_get_cigar(record);

	if ( bam_is_rev(record) ){

		for ( int i = record -> core.n_cigar - 1; i >= 0; i--){

			const int op = bam_cigar_op(cigar[i]); //cigar operation
			const int ol = bam_cigar_oplen(cigar[i]); //number of consecutive operations

			//for a match, advance both reference and query together
			if (op == BAM_CMATCH or op == BAM_CEQUAL or op == BAM_CDIFF){

				for ( int j = refPosition; j < refPosition + ol; j++ ){

					queryPosition++;
				}
				refPosition += ol;
			}
			//for a deletion, advance only the reference position
			else if (op == BAM_CDEL or op == BAM_CREF_SKIP){

				for ( int j = refPosition; j < refPosition + ol; j++ ){

					//std::cout << j << " " << indels.size() << std::endl;
					indels[j] = ol;
				}
				refPosition += ol;
			}
			//for insertions or soft clipping, advance only the query position
			else if (op == BAM_CINS){

				for ( int j = refPosition; j < refPosition + ol; j++ ){

					//std::cout << j << " " << indels.size() << std::endl;
					queryPosition++;
					indels[j] = ol;
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

					queryPosition++;
				}
				refPosition += ol;
			}
			//for a deletion, advance only the reference position
			else if (op == BAM_CDEL or op == BAM_CREF_SKIP){

				for ( int j = refPosition; j < refPosition + ol; j++ ){

					//std::cout << j << " " << indels.size() << std::endl;
					indels[j] = ol;
				}
				refPosition += ol;
			}
			//for insertions or soft clipping, advance only the query position
			else if (op == BAM_CINS){

				for ( int j = refPosition; j < refPosition + ol; j++ ){

					//std::cout << j << " " << indels.size() << std::endl;
					indels[j] = ol;
					queryPosition++;
				}
			}
			//N.B. hard clipping advances neither refernce nor query, so ignore it
		}
	}
	//std::cout << refPosition << " " << indels.size() << std::endl;
	//for (int i = 0; i < indels.size(); i++) std::cout << indels[i] << std::endl;
	return indels;
}


void parseCigar(bam1_t *record, std::map< unsigned int, unsigned int > &ref2query, std::map< unsigned int, unsigned int > &query2ref, int &refStart, int &refEnd ){
	//Covered in: tests/detect/htslib

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
					query2ref[queryPosition] = j;
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
					query2ref[queryPosition] = j;
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
			//N.B. hard clipping advances neither refernce nor query, so ignore it
		}
	}
	refEnd = refStart + refPosition;
}


std::string getQuerySequence( bam1_t *record ){
	//Covered in: tests/detect/htslib

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


void getRefEnd(bam1_t *record, int &refStart, int &refEnd ){
	//Covered in: tests/detect/htslib

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
			//N.B. hard clipping advances neither refernce nor query, so ignore it
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
			//N.B. hard clipping advances neither refernce nor query, so ignore it
		}
	}
	refEnd = refStart + refPosition;
}
