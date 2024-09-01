//----------------------------------------------------------
// Copyright 2019 University of Oxford
// This software is licensed under GPL-3.0.  You should have
// received a copy of the license with this software.  If
// not, please Email the author.
//----------------------------------------------------------

#ifndef READS_H
#define READS_H

#define NFEATURES 5
#define RAWDEPTH 20

#include <string>
#include <map>
#include <vector>
#include <iostream>
#include <cassert>
#include <memory>
#include "htsInterface.h"
#include "common.h"
#include "error_handling.h"
#include "data_IO.h"


struct PoreParameters {

	double shift;
	double scale;
	double eventsPerBase = 0.0;
};


class BandedAlignQCs{

	double avg_log_emission = 0.0;
	bool spanned = false, set = false;
	unsigned int maxGap = 0;
	public:
		void recordQCs(double avg_log_emission, bool spanned, unsigned int maxGap){

			this -> avg_log_emission = avg_log_emission;
			this -> spanned = spanned;
			this -> maxGap = maxGap;
			set = true;
		}
		void printQCs(void){
			assert(set);
			std::cerr << "avg_log_emission" << " " << avg_log_emission << std::endl;
			std::cerr << "spanned" << " " << spanned << std::endl;
			std::cerr << "maxGap" << " " << maxGap << std::endl;
		}
};


struct ReadSegment{
	int leftmostCoord = 0;
	int leftmostIdx = 0;
	int rightmostCoord = 0;
	int rightmostIdx = 0;
	int partners = 0;
	double score = 0.0;
	std::vector<double> stress_signature;
	int querySpan;
};


struct event {

	double mean;
	std::vector< double > raw;
};


class AlignedPosition{

	private:

		bool forTraining = false;
		std::string kmer;
		unsigned int refPos, queryIdx, refIdx;
		std::vector<double> signal;
		double eventAlignQuality;
		std::map<std::string, int> base2index = {{"A",0}, {"T",1}, {"G",2}, {"C",3}};

	public:
		AlignedPosition(std::string kmer, unsigned int refPos, unsigned int queryIdx, unsigned int refIdx, int quality){

			this -> kmer = kmer;
			this -> refPos = refPos;
			this -> queryIdx = queryIdx;
			this -> refIdx = refIdx;
			this -> eventAlignQuality = quality;
		}
		~AlignedPosition() {};
		void addSignal(double s){

			signal.push_back(s);
		}
		std::string getKmer(void){

			return kmer;
		}
		unsigned int getQueryIdx(void){

			return queryIdx;
		}
		unsigned int getReferenceIdx(void){

			return refIdx;
		}
		unsigned int getCoreIndex(void){
			std::string kmer_subseq = kmer.substr(2,5);
			unsigned int kmer_len = kmer_subseq.size();
			unsigned int p = 1;
			unsigned int r = 0;
			for (size_t i = 0; i < kmer_len; i++){

				r += base2index[kmer_subseq.substr(kmer_len-i-1,1)] * p;
				p *= 4;
			}
			r += 1;
			return r;
		}
		unsigned int getResidualIndex(void){
			std::string kmer_subseq = kmer.substr(0,2) + kmer.substr(7,2);
			assert(kmer_subseq.size() == 4);
			unsigned int kmer_len = kmer_subseq.size();
			unsigned int p = 1;
			unsigned int r = 0;
			for (size_t i = 0; i < kmer_len; i++){

				r += base2index[kmer_subseq.substr(kmer_len-i-1,1)] * p;
				p *= 4;
			}
			r += 1;
			return r;
		}
		unsigned int getRefPos(void){

			return refPos;
		}
		double getAlignmentQuality(void){
			
			return eventAlignQuality;
		}
		std::vector<float> makeSignalFeature(void){

			assert(signal.size() > 0);
			assert(kmer.substr(0,1) == "A" || kmer.substr(0,1) == "T" || kmer.substr(0,1) == "G" || kmer.substr(0,1) == "C");
			
			std::vector<float> padded_signal;

			for (size_t i = 0; i < signal.size(); i++){

				padded_signal.push_back(signal[i]);

				if (i == RAWDEPTH - 1) break;
			}
			
			//zero padding if this position has few raw events - these will be masked by the neural network
			if (signal.size() < RAWDEPTH){

				for (size_t i = 0; i < RAWDEPTH - signal.size(); i++){

					padded_signal.push_back(0.);
				}
			}

			assert(padded_signal.size() == RAWDEPTH);
			return padded_signal;
		}
};


namespace DNAscent {

	struct read{
		bam1_t *record;									//alignment record
		std::string basecall;								//basecall sequence (in 5' --> 3' direction on sequencing) from basecaller
		std::string referenceSeqMappedTo;						//subsequence of the reference (in 5' --> 3' direction of sequencing) that the read aligns to
		std::string referenceMappedTo;							//contig name matching the reference
		std::string filename;								//full path to the pod5 or fast5 file containing the raw signal for this read
		std::string humanReadable_detectOut;						//if human readable output is specified, table of analogue calls
		std::string humanReadable_eventalignOut;					//if human readable output is specified, table of aligned events
		std::string readID;								//readID (which may be the result of a split read) from basecaller/MinKNOW
		std::string readID_fetch;							//readID to use for signal fetching from pod5/fast5 - may be equal to readID (if not split) or parent readID (if split)
		PoreParameters scalings;							//shift and scale for signal normalisation
		BandedAlignQCs alignmentQCs;							//quality control measures for the adaptive banded event alignment
		std::vector< event > events;							//downsampled raw signal
		std::vector< double> raw;							//full raw signal in pA
		std::map< unsigned int, unsigned int > refToQuery, queryToRef;			//maps from basecall 0-based indices to referenceSeqMappedTo 0-based indicies (and vice versa)
		std::map< unsigned int, bool > refToDel;					//indicates whether a reference index is in a deletion
		std::vector< std::pair< unsigned int, unsigned int > > eventAlignment;		//rough event alignment from adaptive banded signal alignment
		int refStart, refEnd;								//coordinates on referenceMappedTo (in 5' --> 3' reference direction) where the sequence alignment starts and ends
		bool isReverse = false;								//maps to reverse complement of the reference
		std::string strand = "fwd";							//fwd or rev according to how the read maps
		bool QCpassed = false;								//successfully passed eventalign
		int signalStartCoord = 0;							//if split read: index in the parent signal vector where the signal for this read starts
		int signalLength = -1;								//total length of raw signal as trimmed by Dorado
		int signalTrim = 0;								//number of indices trimmed from the front of the raw signal by Dorado
		std::map<unsigned int, std::shared_ptr<AlignedPosition>> refCoordToAP;		//maps coordinate on the reference contig to an aligned event position
		std::map<unsigned int, std::pair<double,double>> refCoordToCalls;		//maps coordinate on the reference contig to a pair of EdU (first) and BrdU (second) calls
		std::map<unsigned int, std::pair<double,double>> queryIndexToCalls;		//maps index on the query sequence to a pair of EdU (first) and BrdU (second) calls
		size_t pod5_batch;
		size_t pod5_row;
		
		public:
			read(bam1_t *record, bam_hdr_t *bam_hdr, std::map<std::string, IndexEntry> &readID2path, std::map<std::string, std::string> &reference){
				
				this -> record = record;
				
				//get the read name (which will be the ONT readID from basecall)
				const char *queryName = bam_get_qname(record);
				if (queryName == NULL) throw BadBamField("");
				std::string s_queryName(queryName);
				readID = s_queryName;
				readID_fetch = readID; //fetch ID is readID by default unless we have a parent (below)
					
				//signal length
				std::string numSignalsTag = "ns";
				uint8_t *numSignals = bam_aux_get(record, numSignalsTag.c_str());
				signalLength =  bam_aux2i(numSignals);
				
				//coordinates trimmed from signal start
				std::string trimSignalsTag = "ts";
				uint8_t *trimSignal = bam_aux_get(record, trimSignalsTag.c_str());
				signalTrim =  bam_aux2i(trimSignal);
				
				//check if this read was split - if so, get the parentID
				std::string parentIDTag = "pi";
				uint8_t *parentID = bam_aux_get(record, parentIDTag.c_str());
				if (parentID) {
				
					//index in parent signal where the signal for this read starts
					std::string startCoordTag = "sp";
					uint8_t *startCoord = bam_aux_get(record, startCoordTag.c_str());
					signalStartCoord =  bam_aux2i(startCoord);
				
					//get the parent readID and use this one to fetch the raw signal
					char *parentID_char = bam_aux2Z(parentID);
					std::string parentID_str(parentID_char);
					if (parentID_str.size() > 0){
					
						//readID to fetch from pod5 will be the parent readID
						readID_fetch = parentID_str;
					}								
				}
				
				//iterate on the cigar string to fill up the reference-to-query coordinate map
				parseCigar(record, refToQuery, queryToRef, refToDel, refStart, refEnd);

				//get the name of the reference mapped to
				std::string mappedTo(bam_hdr -> target_name[record -> core.tid]);
				referenceMappedTo = mappedTo;

				//unpack index
				IndexEntry ie = readID2path[readID_fetch];
				filename = ie.filepath;
				pod5_batch = ie.batchIndex;
				pod5_row = ie.rowIndex;

				//get the subsequence of the reference this read mapped to
				referenceSeqMappedTo = reference.at(referenceMappedTo).substr(refStart, refEnd - refStart);

				//fetch the basecall from the bam file
				basecall = getQuerySequence(record);

				//account for reverse complements
				if ( bam_is_rev(record) ){

					basecall = reverseComplement( basecall );
					referenceSeqMappedTo = reverseComplement( referenceSeqMappedTo );
					isReverse = true;
					strand = "rev";
				}
			}
			~read(void){ 
				
				bam_destroy1(record); 
			}
			void addSignal(std::string kmer, unsigned int refPos, unsigned int queryIdx, unsigned int refIdx, double sig, int quality){

				if (refCoordToAP.count(refPos) == 0){

					std::shared_ptr<AlignedPosition> ap( new AlignedPosition(kmer, refPos, queryIdx, refIdx, quality));
					ap -> addSignal(sig);
					refCoordToAP[refPos] = ap;
				}
				else{

					refCoordToAP[refPos] -> addSignal(sig);
				}
			}
			std::vector<float> makeSignalTensor(void){

				assert(strand == "fwd" || strand == "rev");
				std::vector<float> tensor;
				tensor.reserve(3 * refCoordToAP.size());

				if (strand == "fwd"){

					for (auto p = refCoordToAP.begin(); p != refCoordToAP.end(); p++){

						std::vector<float> feature = (p -> second) -> makeSignalFeature();
						tensor.insert(tensor.end(), feature.begin(), feature.end());
					}
				}
				else{

					for (auto p = refCoordToAP.rbegin(); p != refCoordToAP.rend(); p++){

						std::vector<float> feature = (p -> second) -> makeSignalFeature();
						tensor.insert(tensor.end(), feature.begin(), feature.end());
					}
				}
				return tensor;
			}
			std::vector<float> makeCoreSequenceTensor(void){

				assert(strand == "fwd" || strand == "rev");
				std::vector<float> tensor;
				tensor.reserve(refCoordToAP.size());

				if (strand == "fwd"){

					for (auto p = refCoordToAP.begin(); p != refCoordToAP.end(); p++){
					
						tensor.push_back( (p -> second) -> getCoreIndex() );
					}
				}
				else{

					for (auto p = refCoordToAP.rbegin(); p != refCoordToAP.rend(); p++){

						tensor.push_back( (p -> second) -> getCoreIndex() );
					}
				}
				return tensor;
			}
			std::vector<float> makeResidualSequenceTensor(void){

				assert(strand == "fwd" || strand == "rev");
				std::vector<float> tensor;
				tensor.reserve(refCoordToAP.size());

				if (strand == "fwd"){

					for (auto p = refCoordToAP.begin(); p != refCoordToAP.end(); p++){
					
						tensor.push_back( (p -> second) -> getResidualIndex() );
					}
				}
				else{

					for (auto p = refCoordToAP.rbegin(); p != refCoordToAP.rend(); p++){

						tensor.push_back( (p -> second) -> getResidualIndex() );
					}
				}
				return tensor;
			}
			std::vector<unsigned int> getReferenceCoords(void){

				std::vector<unsigned int> out;
				out.reserve(refCoordToAP.size());
				if (strand == "fwd"){

					for (auto p = refCoordToAP.begin(); p != refCoordToAP.end(); p++){
						out.push_back(p -> first);
					}
				}
				else{

					for (auto p = refCoordToAP.rbegin(); p != refCoordToAP.rend(); p++){
						out.push_back(p -> first);
					}
				}
				return out;
			}
			std::vector<unsigned int> getQueryIndices(void){

				std::vector<unsigned int> out;
				out.reserve(refCoordToAP.size());
				if (strand == "fwd"){

					for (auto p = refCoordToAP.begin(); p != refCoordToAP.end(); p++){
						out.push_back((p -> second) -> getQueryIdx());
					}
				}
				else{

					for (auto p = refCoordToAP.rbegin(); p != refCoordToAP.rend(); p++){
						out.push_back((p -> second) -> getQueryIdx());
					}
				}
				return out;
			}
			std::vector<unsigned int> getReferenceIndices(void){

				std::vector<unsigned int> out;
				out.reserve(refCoordToAP.size());
				if (strand == "fwd"){

					for (auto p = refCoordToAP.begin(); p != refCoordToAP.end(); p++){
						out.push_back((p -> second) -> getReferenceIdx());
					}
				}
				else{

					for (auto p = refCoordToAP.rbegin(); p != refCoordToAP.rend(); p++){
						out.push_back((p -> second) -> getReferenceIdx());
					}
				}
				return out;
			}
			std::vector<std::string> getKmers(void){

				std::vector<std::string> out;
				out.reserve(refCoordToAP.size());
				if (strand == "fwd"){

					for (auto p = refCoordToAP.begin(); p != refCoordToAP.end(); p++){
						out.push_back((p -> second) -> getKmer());
					}
				}
				else{

					for (auto p = refCoordToAP.rbegin(); p != refCoordToAP.rend(); p++){
						out.push_back((p -> second) -> getKmer());
					}
				}
				return out;
			}
			std::vector<size_t> getSignalShape(void){

				return {refCoordToAP.size(), RAWDEPTH, 1};
			}
			std::vector<size_t> getSequenceShape(void){

				return {refCoordToAP.size()};
			}
			void writeModBamTag(void){

				//get the existing MM tag if there is one
				std::string existing_MMtag = "";
				uint8_t *tagData_MM = bam_aux_get(record, "MM");
				if (tagData_MM != nullptr) {
					existing_MMtag = bam_aux2Z(tagData_MM);
				}
				
				//build the base analogue MM tag for this read and populate vectors of analogue calls
				std::string field_BrdU = "N+b?";
				std::string field_EdU = "N+e?";
				std::vector<uint8_t> BrdUCalls, EdUCalls;
				BrdUCalls.reserve(queryIndexToCalls.size());
				EdUCalls.reserve(queryIndexToCalls.size());				
				unsigned int prevIdx = 0;
				for (auto i = queryIndexToCalls.begin(); i != queryIndexToCalls.end(); i++){
				
					unsigned int queryIndexOfCall = i -> first;
					unsigned int index_delta = queryIndexOfCall - prevIdx;
					
					//append the delta to the string we're building
					field_BrdU += "," + std::to_string(index_delta);
					field_EdU += "," + std::to_string(index_delta);		
					
					prevIdx = queryIndexOfCall + 1;
					
					auto analogueCalls = i -> second;
					EdUCalls.push_back( static_cast<uint8_t>(analogueCalls.first*255.0) );
					BrdUCalls.push_back( static_cast<uint8_t>(analogueCalls.second*255.0) );
				}
				
				//write the MM tag
				std::string tag_MM_value = existing_MMtag + field_BrdU + ";" + field_EdU + ";";
				std::string tagName_MM = "MM";
				bam_aux_append(record, tagName_MM.c_str(), 'Z', int(tag_MM_value.size() + 1), (uint8_t*) tag_MM_value.c_str());				

				//get the existing ML tag if there is one								
				std::vector<uint8_t> probabilities_concat;
				uint8_t *tagData_ML = bam_aux_get(record, "ML");
				if (tagData_ML != nullptr){
				
					uint32_t tagLength_ML = bam_auxB_len(tagData_ML);
					probabilities_concat.reserve(tagLength_ML);
					for (uint32_t i = 0; i < tagLength_ML; i++){
					
						probabilities_concat.push_back( bam_auxB2i(tagData_ML,i) );
					}		
				}
				
				//concatenate the base analogue calls onto the existing tag (if there is one)
				probabilities_concat.insert(probabilities_concat.end(), BrdUCalls.begin(), BrdUCalls.end());
				probabilities_concat.insert(probabilities_concat.end(), EdUCalls.begin(), EdUCalls.end());				

				//write the ML tag
				std::string tagName_ML = "ML";
				bam_aux_update_array(record, tagName_ML.c_str(), 'C', int(probabilities_concat.size()), (uint8_t*) probabilities_concat.data());
			}	
	};
	
	
	struct detectedRead{
		std::string referenceMappedTo;								//contig name matching the reference
		std::string humanReadable_detectOut;							//if human readable output is specified, table of analogue calls
		std::string humanReadable_eventalignOut;						//if human readable output is specified, table of aligned events
		std::string readID;									//readID (which may be the result of a split read) from basecaller/MinKNOW
		std::map< unsigned int, unsigned int > refToQuery, queryToRef;				//maps from basecall 0-based indices to referenceSeqMappedTo 0-based indicies (and vice versa)
		std::map< unsigned int, bool > refToDel;						//indicates whether a reference index is in a deletion
		int refStart, refEnd;									//coordinates on referenceMappedTo (in 5' --> 3' reference direction) where the sequence alignment starts and ends
		bool isReverse = false;									//maps to reverse complement of the reference
		std::string strand = "fwd";								//fwd or rev according to how the read maps
		std::vector< int > referenceCoords;							//coordinates of each call on the reference in 5' --> 3' direction
		std::vector< double > BrdUCalls, EdUCalls;						//probability of BrdU and EdU at each position on referenceCoords (all three vectors same size)
		std::map< int, double > stallScore;							//maps reference coordinate to stall score
		std::vector< int > EdU_segment_label, BrdU_segment_label, thymidine_segment_label;	//binary labels indicating which segment we're in
		std::vector<ReadSegment> EdU_segment, BrdU_segment;					
		std::vector<ReadSegment> origins, terminations, leftForks, rightForks;
		
		public:
			detectedRead(bam1_t *record, bam_hdr_t *bam_hdr){
				
				//get the read name (which will be the ONT readID from basecall)
				const char *queryName = bam_get_qname(record);
				if (queryName == NULL) throw BadBamField("");
				std::string s_queryName(queryName);
				readID = s_queryName;
							
				//iterate on the cigar string to fill up the reference-to-query coordinate map
				parseCigar(record, refToQuery, queryToRef, refToDel, refStart, refEnd);

				//get the name of the reference mapped to
				std::string mappedTo(bam_hdr -> target_name[record -> core.tid]);
				referenceMappedTo = mappedTo;

				//account for reverse complements
				if ( bam_is_rev(record) ){

					isReverse = true;
					strand = "rev";
				}

				//to slice the ML tag
				std::map<std::string, int> fieldToStartIdx, fieldToEndIdx;

				//parse the MM tag and get offsets where each field starts and ends in the ML tag
				uint8_t *tagData_MM = bam_aux_get(record, "MM");
				int offset = 0;
				unsigned int prevQueryIdx = 0;
				if (tagData_MM != nullptr) {
				
					std::string MM_str = bam_aux2Z(tagData_MM);
					std::stringstream ssMM(MM_str);
					std::string analogueField;
					
					//for each ;-delimited analogue field	
					while(std::getline(ssMM, analogueField, ';')){
					
						std::stringstream ssField(analogueField);
						std::string fieldName;
						std::getline(ssField, fieldName, ','); //first one is the field name
						
						std::string keyName;
						if (fieldName == "N+b?") keyName = "BrdU";
						else if (fieldName == "N+e?") keyName = "EdU";
						else keyName = fieldName;
						
						fieldToStartIdx[keyName] = offset;
						
						//for each ,-delimited query skip
						std::string querySkip;
						while(std::getline(ssField, querySkip, ',')){
						
							//WLOG parse the BrdU field to retrieve the query indices
							if (keyName == "BrdU"){
							
								int basesToSkip = stoi(querySkip);
								int indexOnQuery = prevQueryIdx + basesToSkip;
								
								if (queryToRef.count(indexOnQuery) == 0){
								std::cerr << "issue in read" << std::endl;
								}
								
								int indexOnRef = queryToRef.at(indexOnQuery);
								int coordOnRef;
								if (isReverse) coordOnRef = refEnd - indexOnRef;
								else coordOnRef = refStart + indexOnRef;
								
								referenceCoords.push_back(coordOnRef);
								
								prevQueryIdx = indexOnQuery + 1;
							}
							offset++;
						}
						fieldToEndIdx[keyName] = offset;
					}
				}
				
				//fetch the analogue probabilities from the ML tag
				std::vector<double> probabilities_all;
				uint8_t *tagData_ML = bam_aux_get(record, "ML");
				if (tagData_ML != nullptr){
				
					uint32_t tagLength_ML = bam_auxB_len(tagData_ML);
					probabilities_all.reserve(tagLength_ML);
					for (uint32_t i = 0; i < tagLength_ML; i++){
					
						probabilities_all.push_back( static_cast<double>(bam_auxB2i(tagData_ML,i)) / 255.0 );
					}
				}
	
				BrdUCalls = std::vector<double>( probabilities_all.begin() + fieldToStartIdx["BrdU"], probabilities_all.begin() + fieldToEndIdx["BrdU"] );
				EdUCalls = std::vector<double>( probabilities_all.begin() + fieldToStartIdx["EdU"], probabilities_all.begin() + fieldToEndIdx["EdU"] );			

				assert(BrdUCalls.size() == EdUCalls.size());
				assert(BrdUCalls.size() == referenceCoords.size());
				
				if (isReverse){
				
					std::reverse(BrdUCalls.begin(), BrdUCalls.end());
					std::reverse(EdUCalls.begin(), EdUCalls.end());				
					std::reverse(referenceCoords.begin(), referenceCoords.end());								
				}
			}
			detectedRead(std::string readID, std::string contig, int mapLower, int mapUpper, std::string strand, std::vector<int> refCoords, std::vector<double> EdU, std::vector<double> BrdU){
			
				this -> readID = readID.substr(1);
				referenceMappedTo = contig;
				refStart = mapLower;
				refEnd = mapUpper;
				this -> strand = strand;
				if (strand == "rev") isReverse = true;
				referenceCoords = refCoords;
				EdUCalls = EdU;
				BrdUCalls = BrdU;
			}
			std::pair<std::vector<double>, std::vector<double>> getCallFractions(void){

				int resolution = 2000; //look in 2 kb segments
				
				std::vector<double> BrdU_callFractions, EdU_callFractions;
				int BrdUcalls = 0, EdUcalls = 0, attempts = 0, gap = 0, startingPos = -1;
				
				for (size_t i = 0; i < referenceCoords.size(); i++){
			
					if ( BrdUCalls[i] > 0.5 ){
						BrdUcalls++;
						attempts++;
					}
					else if ( EdUCalls[i] > 0.5 ){
						EdUcalls++;
						attempts++;
					}
					else{
						attempts++;
					}

					if ( startingPos == -1 ) startingPos = referenceCoords[i];
					gap = referenceCoords[i] - startingPos;

					if ( gap > resolution and attempts >= resolution / 10 ){

						double BrdU_frac = (double) BrdUcalls / (double) attempts;
						BrdU_callFractions.push_back( BrdU_frac );

						double EdU_frac = (double) EdUcalls / (double) attempts;
						EdU_callFractions.push_back( EdU_frac );

						BrdUcalls = 0, EdUcalls = 0, attempts = 0, gap = 0, startingPos = -1;
					}			
				}
			
				return std::make_pair(BrdU_callFractions, EdU_callFractions);
			}
	};
}


std::vector< std::vector<DNAscent::read *> > sortReadsByFilename(std::vector<DNAscent::read> & );

#endif
