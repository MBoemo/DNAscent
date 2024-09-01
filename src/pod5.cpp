//----------------------------------------------------------
// Copyright 2024 University of Cambridge
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.
// Alternatively, the contents of this file may be used under the terms
// of the GNU General Public License v3 or later.
//----------------------------------------------------------


#include "pod5.h"
#include "../pod5-file-format/c++/pod5_format/c_api.h"
#include "error_handling.h"
#include <vector>
#include <string>
#include <iostream>
#include <array>
#include <boost/lexical_cast.hpp>
#include <boost/uuid/uuid.hpp>
#include <boost/uuid/uuid_io.hpp>


//adapted from https://github.com/nanoporetech/pod5-file-format
void pod5_getSignal( DNAscent::read &r ){
	
	std::string filepath = r.filename;
	size_t batch_index = r.pod5_batch;
	size_t batch_row = r.pod5_row;
	
	Pod5FileReader_t *file = pod5_open_file(filepath.c_str());
	if (!file) {
		std::cerr << "Failed to open file " << filepath << ": " << pod5_get_error_string() << "\n";
		std::cerr << "readID:  " << r.readID_fetch << "\n";
		throw BadPod5Field();
	}

	Pod5ReadRecordBatch_t *batch = nullptr;
	if (pod5_get_read_batch(&batch, file, batch_index) != POD5_OK) {
		std::cerr << "Failed to get batch: " << pod5_get_error_string() << "\n";
		throw BadPod5Field();
	}

	uint16_t read_table_version = 0;
	ReadBatchRowInfo_t read_data;
	if (pod5_get_read_batch_row_info_data(batch, batch_row, READ_BATCH_ROW_INFO_VERSION, &read_data, &read_table_version) != POD5_OK){
		std::cerr << "Failed to get read " << batch_row << "\n";
		throw BadPod5Field();
	}

	std::size_t sample_count = 0;
	pod5_get_read_complete_sample_count(file, batch, batch_row, &sample_count);

	std::vector<std::int16_t> samples;
	samples.resize(sample_count);
	pod5_get_read_complete_signal(file, batch, batch_row, samples.size(), samples.data());
	
	//normalise signal to pA
	r.raw.reserve(sample_count);
	for ( size_t i = 0; i < sample_count; i++ ){
		r.raw.push_back( ((float) samples[i] + (float) read_data.calibration_offset) * (float) read_data.calibration_scale );
	}
	
	//trim raw signal from the parent if this is a split read
	if (r.readID != r.readID_fetch){
	
		assert(r.signalLength > 0);
		size_t sig_start = r.signalStartCoord + r.signalTrim;
		size_t sig_end = r.signalStartCoord + r.signalLength;
		r.raw.erase(r.raw.begin(), r.raw.begin() + sig_start);
		r.raw.erase(r.raw.begin() + (sig_end - sig_start), r.raw.end());
	}
	else{ //trim using the normal dorado bounds if the read wasn't split

		size_t sig_start = r.signalTrim;
		size_t sig_end = r.signalLength;
		r.raw.erase(r.raw.begin(), r.raw.begin() + sig_start);
		r.raw.erase(r.raw.begin() + (sig_end - sig_start), r.raw.end());
	}

	if (pod5_free_read_batch(batch) != POD5_OK) {
		std::cerr << "Failed to release batch\n";
		throw BadPod5Field();
	}
	
	// Close the reader
	if (pod5_close_and_free_reader(file) != POD5_OK) {
		std::cerr << "Failed to close reader: " << pod5_get_error_string() << "\n";
		throw BadPod5Field();
	}
}


//adapted from https://github.com/nanoporetech/pod5-file-format
void pod5_getSignal_batch( std::vector<DNAscent::read *> readBatch ){
	
	std::string filepath = readBatch[0] -> filename;
	std::vector<boost::uuids::uuid> search_uuids;
	std::map<std::string, std::vector<DNAscent::read *>> fetchIDToRead;
	
	for (size_t i = 0; i < readBatch.size(); i++){
	
		search_uuids.push_back(boost::lexical_cast<boost::uuids::uuid>(readBatch[i] -> readID_fetch));
		
		if (fetchIDToRead.count(readBatch[i] -> readID_fetch) > 0){
		
			fetchIDToRead[readBatch[i] -> readID_fetch].push_back(readBatch[i]);
		}
		else{
			fetchIDToRead[readBatch[i] -> readID_fetch] ={readBatch[i]};		
		}
		
		//all reads should be from the same pod5 file
		assert(readBatch[i] -> filename == filepath);
	}
	
	pod5_init();
	
	Pod5FileReader_t *file = pod5_open_file(filepath.c_str());

	std::size_t batch_count = 0;
	if (pod5_get_read_batch_count(&batch_count, file) != POD5_OK) {
		std::cerr << "Failed to query batch count: " << pod5_get_error_string() << "\n";
		throw BadPod5Field();
	}
	
	std::vector<std::uint32_t> traversal_batch_counts(batch_count);
	std::vector<std::uint32_t> traversal_row_indices(search_uuids.size());
	std::size_t find_success_count = 0;
	if (pod5_plan_traversal(file, (uint8_t *)search_uuids.data(), search_uuids.size(), traversal_batch_counts.data(), traversal_row_indices.data(), &find_success_count) != POD5_OK){
		std::cerr << "Failed to plan traversal of file: " << pod5_get_error_string() << "\n";
		throw BadPod5Field();
	}
	
	if (find_success_count != search_uuids.size()) {
		std::cerr << "Failed to find " << (search_uuids.size() - find_success_count) << " reads\n";
		throw BadPod5Field();
	}

	std::size_t read_count = 0;
	std::size_t row_offset = 0;	

	for (size_t batch_index = 0; batch_index < batch_count; batch_index++) {
		Pod5ReadRecordBatch_t *batch = nullptr;
		if (pod5_get_read_batch(&batch, file, batch_index) != POD5_OK) {
			std::cerr << "Failed to get batch: " << pod5_get_error_string() << "\n";
			throw BadPod5Field();
		}

		std::size_t batch_row_count = 0;
		if (pod5_get_read_batch_row_count(&batch_row_count, batch) != POD5_OK) {
			std::cerr << "Failed to get batch row count\n";
			throw BadPod5Field();
		}
		
		for (std::size_t row_index = 0; row_index < traversal_batch_counts[batch_index];++row_index){
			std::uint32_t batch_row = traversal_row_indices[row_index + row_offset];

			uint16_t read_table_version = 0;
			ReadBatchRowInfo_t read_data;
			if (pod5_get_read_batch_row_info_data(batch, batch_row, READ_BATCH_ROW_INFO_VERSION, &read_data, &read_table_version) != POD5_OK){
				std::cerr << "Failed to get read " << batch_row << "\n";
				throw BadPod5Field();
			}

			std::size_t sample_count = 0;
			pod5_get_read_complete_sample_count(file, batch, batch_row, &sample_count);

			std::vector<std::int16_t> samples;
			samples.resize(sample_count);
			pod5_get_read_complete_signal(file, batch, batch_row, samples.size(), samples.data());
			
			//normalise signal to pA
			std::vector<double> signal_pA;
			signal_pA.reserve(sample_count);
			for ( size_t i = 0; i < sample_count; i++ ){
				signal_pA.push_back( ((float) samples[i] + (float) read_data.calibration_offset) * (float) read_data.calibration_scale );
			}
			
			//add the signal to the read (or reads) in this batch that match the fetchID
			std::array<char, 37> formatted_read_id;
			pod5_format_read_id(read_data.read_id, formatted_read_id.data());
			std::string thisID = formatted_read_id.data();
			
			for (auto r = fetchIDToRead[thisID].begin(); r < fetchIDToRead[thisID].end(); r++){
			
				//trim raw signal from the parent if this is a split read
				if ((*r) -> readID != (*r) -> readID_fetch){
					assert((*r) -> signalLength > 0);
					signal_pA = std::vector<double>(signal_pA.begin() + (*r) -> signalStartCoord + (*r) -> signalTrim,
					                                signal_pA.begin() + (*r) -> signalStartCoord + (*r) -> signalLength);
				}
				else{ //trim using the normal dorado bounds if the read wasn't split
					signal_pA = std::vector<double>(signal_pA.begin() + (*r) -> signalTrim,
					                                signal_pA.begin() + (*r) -> signalLength);				
				}
			
				(*r) -> raw = signal_pA;
			}
			
			read_count++;
		}
		row_offset += traversal_batch_counts[batch_index];

		if (pod5_free_read_batch(batch) != POD5_OK) {
			std::cerr << "Failed to release batch\n";
			throw BadPod5Field();
		}
	}

	// Close the reader
	if (pod5_close_and_free_reader(file) != POD5_OK) {
		std::cerr << "Failed to close reader: " << pod5_get_error_string() << "\n";
		throw BadPod5Field();
	}

	// Cleanup the library
	pod5_terminate();
}


//adapted from https://github.com/nanoporetech/pod5-file-format
std::vector< std::string > pod5_extract_readIDs(std::string filepath){

	std::vector< std::string > readIDs;
	
	pod5_init();
	
	Pod5FileReader_t *file = pod5_open_file(filepath.c_str());

	std::size_t batch_count = 0;
	if (pod5_get_read_batch_count(&batch_count, file) != POD5_OK) {
		std::cerr << "Failed to query batch count: " << pod5_get_error_string() << "\n";
		throw BadPod5Field();
	}

	std::size_t read_count = 0;

	for (size_t batch_index = 0; batch_index < batch_count; batch_index++) {
		Pod5ReadRecordBatch_t *batch = nullptr;
		if (pod5_get_read_batch(&batch, file, batch_index) != POD5_OK) {
			std::cerr << "Failed to get batch: " << pod5_get_error_string() << "\n";
			throw BadPod5Field();
		}

		std::size_t batch_row_count = 0;
		if (pod5_get_read_batch_row_count(&batch_row_count, batch) != POD5_OK) {
			std::cerr << "Failed to get batch row count\n";
			throw BadPod5Field();
		}

		for (std::size_t row = 0; row < batch_row_count; row++) {
			uint16_t read_table_version = 0;
			ReadBatchRowInfo_t read_data;
			if (pod5_get_read_batch_row_info_data(batch, row, READ_BATCH_ROW_INFO_VERSION, &read_data, &read_table_version) != POD5_OK){
				std::cerr << "Failed to get read " << row << "\n";
				throw BadPod5Field();
			}

			std::array<char, 37> formatted_read_id;
			pod5_format_read_id(read_data.read_id, formatted_read_id.data());
			
			std::string thisID = formatted_read_id.data();
			thisID += "\t" + std::to_string(batch_index) + "\t" + std::to_string(row);
			
			readIDs.push_back(thisID);
			read_count += 1;
		}

		if (pod5_free_read_batch(batch) != POD5_OK) {
			std::cerr << "Failed to release batch\n";
			throw BadPod5Field();
		}
	}

	// Close the reader
	if (pod5_close_and_free_reader(file) != POD5_OK) {
		std::cerr << "Failed to close reader: " << pod5_get_error_string() << "\n";
		throw BadPod5Field();
	}

	// Cleanup the library
	pod5_terminate();
	
	return readIDs;

}
