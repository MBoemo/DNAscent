//----------------------------------------------------------
// Copyright 2019-2020 University of Oxford
// This software is licensed under GPL-3.0.  You should have
// received a copy of the license with this software.  If
// not, please Email the author.
//----------------------------------------------------------


#include "../fast5/include/fast5.hpp"
#include "scrappie/event_detection.h"
#include "scrappie/scrappie_common.h"
#include "error_handling.h"
#include "fast5.h"
#include <vector>
#include <string>

#define _USE_MATH_DEFINES

//from scrappie
float fast5_read_float_attribute(hid_t group, const char *attribute) {
	float val = NAN;
	if (group < 0) {
#ifdef DEBUG_FAST5_IO
		fprintf(stderr, "Invalid group passed to %s:%d.", __FILE__, __LINE__);
#endif
		return val;
	}

	hid_t attr = H5Aopen(group, attribute, H5P_DEFAULT);
	if (attr < 0) {
#ifdef DEBUG_FAST5_IO
		fprintf(stderr, "Failed to open attribute '%s' for reading.", attribute);
#endif
		return val;
	}

	H5Aread(attr, H5T_NATIVE_FLOAT, &val);
	H5Aclose(attr);

	return val;
}
//end scrappie


void fast5_getSignal( DNAscent::read &r ){

	std::string filename = r.filename;
	std::string readID = r.readID_fetch;

	//open the file
	hid_t hdf5_file = H5Fopen(filename.c_str(), H5F_ACC_RDONLY, H5P_DEFAULT);
	if (hdf5_file < 0) throw IOerror(filename.c_str());

	//check for vbz compression and fail if plugin is not loaded
	bool vbz_compressed = false;
    unsigned int flags;
	size_t nelmts = 1;
    unsigned int values_out[1] = {99};
	char filter_name[80];
	std::string signal_path = "/read_" + readID + "/Raw/Signal";
	hid_t dcheck = H5Dopen(hdf5_file, signal_path.c_str(), H5P_DEFAULT);
    hid_t dcpl = H5Dget_create_plist(dcheck);
    H5Z_filter_t filter_id = H5Pget_filter2(dcpl, (unsigned) 0, &flags, &nelmts, values_out, sizeof(filter_name) - 1, filter_name, NULL);
    H5Pclose (dcpl);
    H5Dclose (dcheck);
	if(filter_id == 32020) vbz_compressed = true;

	//get the channel parameters
	std::string scaling_path = "/read_" + readID + "/channel_id";
	hid_t scaling_group = H5Gopen(hdf5_file, scaling_path.c_str(), H5P_DEFAULT);
	float digitisation = fast5_read_float_attribute(scaling_group, "digitisation");
	float offset = fast5_read_float_attribute(scaling_group, "offset");
	float range = fast5_read_float_attribute(scaling_group, "range");
	H5Gclose(scaling_group);

	//get the raw signal
	hid_t space;
	hsize_t nsample;
	float raw_unit;
	float *rawptr = NULL;

	hid_t dset = H5Dopen(hdf5_file, signal_path.c_str(), H5P_DEFAULT);
	if (dset < 0 ) throw BadFast5Field(); 
	space = H5Dget_space(dset);
	if (space < 0 ) throw BadFast5Field(); 
	H5Sget_simple_extent_dims(space, &nsample, NULL);
   	rawptr = (float*)calloc(nsample, sizeof(float));
    herr_t status = H5Dread(dset, H5T_NATIVE_FLOAT, H5S_ALL, H5S_ALL, H5P_DEFAULT, rawptr);

	if ( status < 0 ){
		free(rawptr);
		H5Dclose(dset);

		if (vbz_compressed) throw VBZError(filename);

		throw BadFast5Field();
	}
	H5Dclose(dset);

	std::vector<double> signal_pA;
	signal_pA.reserve(nsample);
	
	raw_unit = range / digitisation;
	for ( size_t i = 0; i < nsample; i++ ){

		signal_pA.push_back( (rawptr[i] + offset) * raw_unit );
	}
	
	r.raw = signal_pA;

	//fail on empty signal
	if (r.raw.size() == 0){

		std::cerr << "Empty signal found in fast5 file." << std::endl;
		std::cerr << "   ReadID: " << r.readID << std::endl;
		std::cerr << "   Filename: " << r.filename << std::endl;
		exit(EXIT_FAILURE);
	}

	free(rawptr);
	
	H5Fclose(hdf5_file);
}


void fast5_getSignal_batch( std::vector<DNAscent::read *> readBatch ){

	std::string filename = readBatch[0] -> filename;

	//open the file
	hid_t hdf5_file = H5Fopen(filename.c_str(), H5F_ACC_RDONLY, H5P_DEFAULT);
	if (hdf5_file < 0) throw IOerror(filename.c_str());
	
	for (size_t i = 0; i < readBatch.size(); i++){
	
		std::string readID = readBatch[i] -> readID_fetch;

		//get the channel parameters
		std::string scaling_path = "/read_" + readID + "/channel_id";
		hid_t scaling_group = H5Gopen(hdf5_file, scaling_path.c_str(), H5P_DEFAULT);
		float digitisation = fast5_read_float_attribute(scaling_group, "digitisation");
		float offset = fast5_read_float_attribute(scaling_group, "offset");
		float range = fast5_read_float_attribute(scaling_group, "range");
		H5Gclose(scaling_group);

		//get the raw signal
		hid_t space;
		hsize_t nsample;
		float raw_unit;
		float *rawptr = NULL;

		std::string signal_path = "/read_" + readID + "/Raw/Signal";
		hid_t dset = H5Dopen(hdf5_file, signal_path.c_str(), H5P_DEFAULT);
		if (dset < 0 ) throw BadFast5Field(); 
		space = H5Dget_space(dset);
		if (space < 0 ) throw BadFast5Field(); 
		H5Sget_simple_extent_dims(space, &nsample, NULL);
	   	rawptr = (float*)calloc(nsample, sizeof(float));
	    	herr_t status = H5Dread(dset, H5T_NATIVE_FLOAT, H5S_ALL, H5S_ALL, H5P_DEFAULT, rawptr);

		if ( status < 0 ){
			free(rawptr);
			H5Dclose(dset);
			throw BadFast5Field();
		}
		H5Dclose(dset);

		std::vector<double> signal_pA;
		signal_pA.reserve(nsample);
		
		raw_unit = range / digitisation;
		for ( size_t i = 0; i < nsample; i++ ){

			signal_pA.push_back( (rawptr[i] + offset) * raw_unit );
		}
		
		readBatch[i] -> raw = signal_pA;

		free(rawptr);
	}
	
	H5Fclose(hdf5_file);
}

std::vector<std::string> fast5_extract_readIDs(std::string filepath){

    hid_t hdf5_file = H5Fopen(filepath.c_str(), H5F_ACC_RDONLY, H5P_DEFAULT);
    if (hdf5_file < 0) throw IOerror(filepath.c_str());

    std::vector<std::string> out;
    ssize_t buffer_size = 0;
    char* buffer = NULL;

    // get the number of groups in the root group
    H5G_info_t group_info;
    int ret = H5Gget_info_by_name(hdf5_file, "/", &group_info, H5P_DEFAULT);
    if(ret < 0) {
        fprintf(stderr, "error getting group info\n");
        exit(EXIT_FAILURE);
    }

    for(size_t group_idx = 0; group_idx < group_info.nlinks; ++group_idx) {

        // retrieve the size of this group name
        ssize_t size = H5Lget_name_by_idx(hdf5_file, "/", H5_INDEX_NAME, H5_ITER_INC, group_idx, NULL, 0, H5P_DEFAULT);
	
        if(size < 0) {
            fprintf(stderr, "error getting group name size\n");
            exit(EXIT_FAILURE);
        }
        size += 1; // for null terminator
           
        if(size > buffer_size) {
            buffer = (char*)realloc(buffer, size);
            buffer_size = size;
        }
    
        // copy the group name
        H5Lget_name_by_idx(hdf5_file, "/", H5_INDEX_NAME, H5_ITER_INC, group_idx, buffer, buffer_size, H5P_DEFAULT);
        buffer[size] = '\0';
        
        char prefix[] = "read_";
	size_t len = strlen(prefix);
	if (strncmp(buffer, prefix, len) == 0) {
		memmove(buffer, buffer + len, strlen(buffer + len) + 1);
	}       
        
        out.push_back(buffer);
	
    }

    free(buffer);
    buffer = NULL;
    buffer_size = 0;
    return out;
}
