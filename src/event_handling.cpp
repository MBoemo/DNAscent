//----------------------------------------------------------
// Copyright 2019-2020 University of Oxford
// Written by Michael A. Boemo (mb915@cam.ac.uk)
// This software is licensed under GPL-3.0.  You should have
// received a copy of the license with this software.  If
// not, please Email the author.
//----------------------------------------------------------

//#define EVENT_LENGTHS 1
//#define SHOW_PROGRESS 1

#include <iterator>
#include <algorithm>
#include <math.h>
#include "probability.h"
#include "error_handling.h"
#include "event_handling.h"
#include "../fast5/include/fast5.hpp"
#include <chrono>
#include "config.h"

#include "scrappie/event_detection.h"
#include "scrappie/scrappie_common.h"

// #include "../pod5-file-format/c++/pod5_format/c_api.h"


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


void fast5_getSignal( read &r ){

	//open the file
	hid_t hdf5_file = H5Fopen(r.filename.c_str(), H5F_ACC_RDONLY, H5P_DEFAULT);
	if (hdf5_file < 0) throw IOerror(r.filename.c_str());

	//get the channel parameters
	std::string scaling_path = "/read_" + r.readID + "/channel_id";
	hid_t scaling_group = H5Gopen(hdf5_file, scaling_path.c_str(), H5P_DEFAULT);
	float digitisation = fast5_read_float_attribute(scaling_group, "digitisation");
	float offset = fast5_read_float_attribute(scaling_group, "offset");
	float range = fast5_read_float_attribute(scaling_group, "range");
	//float sample_rate = fast5_read_float_attribute(scaling_group, "sampling_rate");
	H5Gclose(scaling_group);

	//get the raw signal
	hid_t space;
	hsize_t nsample;
	float raw_unit;
	float *rawptr = NULL;

	std::string signal_path = "/read_" + r.readID + "/Raw/Signal";
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
		return;
	}
	H5Dclose(dset);
	
	raw_unit = range / digitisation;
	r.raw.reserve(nsample);
	for ( size_t i = 0; i < nsample; i++ ){

		r.raw.push_back( (rawptr[i] + offset) * raw_unit );
	}

	free(rawptr);
	H5Fclose(hdf5_file);
}


PoreParameters estimateScaling_theilSen(std::vector< double > &signals, std::vector<unsigned int> &kmer_ranks, PoreParameters s, bool useFitPoreModel){

	assert(signals.size() == kmer_ranks.size());

	size_t maxPoints = 1000;
	size_t trimSize = 50;
	size_t minLength = maxPoints;
	
	//for short reads, exit without doing refinement of scaling parameters	
	if (kmer_ranks.size() < minLength) return s;

	size_t effectiveSize = signals.size() - 2*trimSize;
	
	size_t skipInterval = 1;
	size_t numPoints = effectiveSize;
	if (effectiveSize > maxPoints){
		skipInterval = effectiveSize/maxPoints;
		numPoints = maxPoints;
	}
	
	std::vector<double> x, y;

	size_t i = trimSize;
	for (size_t j = 0; j < numPoints; j++){
	
		assert(i < signals.size() and i < kmer_ranks.size());
	
		x.push_back( (signals[i]-s.shift)/s.scale );

		std::pair<double,double> meanStd;
		if (useFitPoreModel){
			meanStd = Pore_Substrate_Config.unlabelled_model[kmer_ranks[i]];		
		}
		else{
			meanStd = Pore_Substrate_Config.pore_model[kmer_ranks[i]];
		}	
		double kmer_mean = meanStd.first;
		y.push_back(kmer_mean);
		i += skipInterval;
	}

	std::vector<double> slopes;
	slopes.reserve(maxPoints * maxPoints);
	for (size_t i = 0; i < x.size(); i++){	
		for (size_t j = i+1; j < x.size(); j++){
			
			double dy = y[i] - y[j];
			double dx = x[i] - x[j];

			slopes.push_back( dy /dx );
		}	
	}
	
	std::sort(slopes.begin(), slopes.end());
	double slope_median = slopes[ slopes.size() / 2 ];
	std::vector<double> intercepts;
	intercepts.reserve(x.size());
	for (size_t i = 0; i < x.size(); i++){	

		intercepts.push_back( y[i] - slope_median*x[i] );
	}

	std::sort(intercepts.begin(), intercepts.end());
	double intercept_median = intercepts[ intercepts.size() / 2 ];
	PoreParameters params_ts;
	
	if (slope_median == 0.){
	
		params_ts.shift = -1.;
		params_ts.scale = -1.;
		return params_ts;		
	}
	
	//use TS parameters to refine shift and scale similar to Remora
	double scale_correlation = 1. / slope_median;
	double shift_correlation = -intercept_median / slope_median;
	double shift_TSrefined = s.shift + (shift_correlation * s.scale);
	double scale_TSrefined = s.scale * scale_correlation;

	//std::cerr << "ROUGH TO TS: " << s.scale << " " << scale_TSrefined << " " << s.shift << " " << shift_TSrefined << std::endl;
	//std::cerr << s.scale << " " << scale_TSrefined << " " << s.shift << " " << shift_TSrefined << std::endl;

	params_ts.shift = shift_TSrefined;
	params_ts.scale = scale_TSrefined;

	return params_ts;
}


//start: adapted from nanopolish (https://github.com/jts/nanopolish)
//licensed under MIT

inline float logProbabilityMatch(unsigned int kmerIndex, event e, double shift, double scale, bool useFitPoreModel){

	std::pair<double,double> meanStd;
	if (useFitPoreModel){
		meanStd = Pore_Substrate_Config.unlabelled_model[kmerIndex];		
	}
	else{
		meanStd = Pore_Substrate_Config.pore_model[kmerIndex];
	}

	double mu = meanStd.first;
	double sigma = meanStd.second;
	
	//scale the signal to the pore model
	double x = (e.mean - shift)/scale;
	
	//normal distribution
	float a = (x - mu) / sigma;	
	static const float log_inv_sqrt_2pi = log(0.3989422804014327);
	double thymProb = log_inv_sqrt_2pi - eln(sigma) + (-0.5f * a * a);
	return thymProb;
}	

#define event_kmer_to_band(ei, ki) (ei + 1) + (ki + 1)
#define band_event_to_offset(bi, ei) band_lower_left[bi].event_idx - (ei)
#define band_kmer_to_offset(bi, ki) (ki) - band_lower_left[bi].kmer_idx
#define is_offset_valid(offset) (offset) >= 0 && (offset) < bandwidth
#define event_at_offset(bi, offset) band_lower_left[(bi)].event_idx - (offset)
#define kmer_at_offset(bi, offset) band_lower_left[(bi)].kmer_idx + (offset)
#define move_down(curr_band) { curr_band.event_idx + 1, curr_band.kmer_idx }
#define move_right(curr_band) { curr_band.event_idx, curr_band.kmer_idx + 1 }

std::pair<std::vector<double>, std::vector<unsigned int>> adaptive_banded_simple_event_align( read &r, std::vector<unsigned int> &kmer_ranks_query, std::vector<unsigned int> &kmer_ranks_ref, bool useFitPoreModel ){

	//benchmarking
	//std::chrono::steady_clock::time_point tp1 = std::chrono::steady_clock::now();

	std::string sequence = r.basecall;

	size_t k = Pore_Substrate_Config.kmer_len;
	size_t n_events = r.events.size();
	size_t n_kmers = sequence.size() - k + 1;

	// backtrack markers
	const uint8_t FROM_D = 0;
	const uint8_t FROM_U = 1;
	const uint8_t FROM_L = 2;
 
	// qc
	double min_average_log_emission = Pore_Substrate_Config.AdaptiveBanded_config.min_average_log_emission;
	int max_gap_threshold = Pore_Substrate_Config.AdaptiveBanded_config.max_gap_threshold;

	// banding
	int bandwidth = Pore_Substrate_Config.AdaptiveBanded_config.bandwidth;

	int half_bandwidth = bandwidth / 2;
 
	// transition penalties
	double events_per_kmer = (double)n_events / n_kmers;
	double p_stay = 1 - (1 / (events_per_kmer + 1));

	// setting a tiny skip penalty helps keep the true alignment within the adaptive band
	// this was empirically determined
	double epsilon = 1e-30;
	double lp_skip = log(epsilon);
	double lp_stay = log(p_stay);
	double lp_step = log(1.0 - exp(lp_skip) - exp(lp_stay));
	double lp_trim = log(0.01);
 
	// dp matrix
	size_t n_rows = n_events + 1;
	size_t n_cols = n_kmers + 1;
	size_t n_bands = n_rows + n_cols;
 
	// Initialize
	typedef std::vector<float> bandscore;
	typedef std::vector<uint8_t> bandtrace;

	std::vector<bandscore> bands(n_bands);
	std::vector<bandtrace> trace(n_bands);

	for(size_t i = 0; i < n_bands; ++i) {
		bands[i].resize(bandwidth, -INFINITY);
		trace[i].resize(bandwidth, 0);
	}

	// Keep track of the event/kmer index for the lower left corner of the band
	// these indices are updated at every iteration to perform the adaptive banding
	// Only the first two bands have their coordinates initialized, the rest are computed adaptively
	struct EventKmerPair {
		int event_idx;
		int kmer_idx;
	};

	std::vector<EventKmerPair> band_lower_left(n_bands);
 
	// initialize range of first two bands
	band_lower_left[0].event_idx = half_bandwidth - 1;
	band_lower_left[0].kmer_idx = -1 - half_bandwidth;
	band_lower_left[1] = move_down(band_lower_left[0]);

	// band 0: score zero in the central cell
	int start_cell_offset = band_kmer_to_offset(0, -1);
	assert(is_offset_valid(start_cell_offset));
	assert(band_event_to_offset(0, -1) == start_cell_offset);
	bands[0][start_cell_offset] = 0.0f;
    
	// band 1: first event is trimmed
	int first_trim_offset = band_event_to_offset(1, 0);
	assert(kmer_at_offset(1, first_trim_offset) == -1);
	assert(is_offset_valid(first_trim_offset));
	bands[1][first_trim_offset] = lp_trim;
	trace[1][first_trim_offset] = FROM_U;

	int fills = 0;

	// fill in remaining bands
	for(unsigned int band_idx = 2; band_idx < n_bands; ++band_idx) {
	// Determine placement of this band according to Suzuki's adaptive algorithm
        // When both ll and ur are out-of-band (ob) we alternate movements
        // otherwise we decide based on scores
		float ll = bands[band_idx - 1][0];
		float ur = bands[band_idx - 1][bandwidth - 1];
		bool ll_ob = ll == -INFINITY;
		bool ur_ob = ur == -INFINITY;
        
		bool right = false;
		if(ll_ob && ur_ob) {
			right = band_idx % 2 == 1;
		} else {
			right = ll < ur; // Suzuki's rule
		}

		if(right) {
			band_lower_left[band_idx] = move_right(band_lower_left[band_idx - 1]);
		} else {
			band_lower_left[band_idx] = move_down(band_lower_left[band_idx - 1]);
		}

		// If the trim state is within the band, fill it in here
		int trim_offset = band_kmer_to_offset(band_idx, -1);
		if(is_offset_valid(trim_offset)) {
			unsigned int event_idx = event_at_offset(band_idx, trim_offset);
			if(event_idx >= 0 && event_idx < n_events) {
				bands[band_idx][trim_offset] = lp_trim * (event_idx + 1);
				trace[band_idx][trim_offset] = FROM_U;
			} else {
				bands[band_idx][trim_offset] = -INFINITY;
			}
		}

		// Get the offsets for the first and last event and kmer
		// We restrict the inner loop to only these values
		int kmer_min_offset = band_kmer_to_offset(band_idx, 0);
		int kmer_max_offset = band_kmer_to_offset(band_idx, n_kmers);
		int event_min_offset = band_event_to_offset(band_idx, n_events - 1);
		int event_max_offset = band_event_to_offset(band_idx, -1);

		int min_offset = std::max(kmer_min_offset, event_min_offset);
		min_offset = std::max(min_offset, 0);

		int max_offset = std::min(kmer_max_offset, event_max_offset);
		max_offset = std::min(max_offset, bandwidth);

		for(int offset = min_offset; offset < max_offset; ++offset) {
			int event_idx = event_at_offset(band_idx, offset);
			int kmer_idx = kmer_at_offset(band_idx, offset);

			unsigned int kmer_rank = kmer_ranks_query[kmer_idx];
 
			int offset_up   = band_event_to_offset(band_idx - 1, event_idx - 1); 
			int offset_left = band_kmer_to_offset(band_idx - 1, kmer_idx - 1);
			int offset_diag = band_kmer_to_offset(band_idx - 2, kmer_idx - 1);

			float up   = is_offset_valid(offset_up)   ? bands[band_idx - 1][offset_up]   : -INFINITY;
			float left = is_offset_valid(offset_left) ? bands[band_idx - 1][offset_left] : -INFINITY;
			float diag = is_offset_valid(offset_diag) ? bands[band_idx - 2][offset_diag] : -INFINITY;
 
			float lp_emission = logProbabilityMatch(kmer_rank, r.events[event_idx], r.scalings.shift, r.scalings.scale, useFitPoreModel);

			float score_d = diag + lp_step + lp_emission;
			float score_u = up + lp_stay + lp_emission;
			float score_l = left + lp_skip;

			float max_score = score_d;
			uint8_t from = FROM_D;

			max_score = score_u > max_score ? score_u : max_score;
			from = max_score == score_u ? FROM_U : from;
			max_score = score_l > max_score ? score_l : max_score;
			from = max_score == score_l ? FROM_L : from;
	    
			bands[band_idx][offset] = max_score;
			trace[band_idx][offset] = from;
			fills += 1;
		}
	}

	//benchmarking
	//std::chrono::steady_clock::time_point tp2 = std::chrono::steady_clock::now();
	//std::cout << "banded alignment: " << std::chrono::duration_cast<std::chrono::microseconds>(tp2 - tp1).count() << std::endl;

	//
	// Backtrack to compute alignment
	//
	double sum_emission = 0.;
	double n_aligned_events = 0;
    
	float max_score = -INFINITY;
	int curr_event_idx = 0;
	int curr_kmer_idx = n_kmers -1;

	// Find best score between an event and the last k-mer. after trimming the remaining evnets
	for(unsigned int event_idx = 0; event_idx < n_events; ++event_idx) {
		unsigned int band_idx = event_kmer_to_band(event_idx, curr_kmer_idx);
		assert(band_idx < bands.size());
		int offset = band_event_to_offset(band_idx, event_idx);
		if(is_offset_valid(offset)) {
			float s = bands[band_idx][offset] + (n_events - event_idx) * lp_trim;
			if(s > max_score) {
				max_score = s;
				curr_event_idx = event_idx;
			}
		}
	}
//end adapted from nanopolish

	//benchmarking
	//std::chrono::steady_clock::time_point tp3 = std::chrono::steady_clock::now();
	//std::cout << "calculate end index: " << std::chrono::duration_cast<std::chrono::microseconds>(tp3 - tp2).count() << std::endl;

	r.eventAlignment.reserve(r.events.size());

	int curr_gap = 0;
	int max_gap = 0;

	std::vector< double > signalBuffer;
	std::vector< unsigned int > cleanedRanks;
	std::vector< double > cleanedSignals;

	while(curr_kmer_idx >= 0 && curr_event_idx >= 0) {
        
		// emit alignment
		r.eventAlignment.push_back(std::make_pair(curr_event_idx, curr_kmer_idx));

		// qc stats
		unsigned int kmer_rank = kmer_ranks_query[curr_kmer_idx];
		float logProbability = logProbabilityMatch(kmer_rank, r.events[curr_event_idx], r.scalings.shift, r.scalings.scale, useFitPoreModel);
		sum_emission += logProbability;

		n_aligned_events += 1;
		
		//TESTING - print the alignment
		//double mu = meanStd.first;
		//double sigma = meanStd.second;
		//double x = (r.events[curr_event_idx].mean - r.scalings.shift)/r.scalings.scale;
		//std::cout << curr_kmer_idx << "\t" << x << "\t" << mu << std::endl;
		//ENDTESTING

		int band_idx = event_kmer_to_band(curr_event_idx, curr_kmer_idx);
		int offset = band_event_to_offset(band_idx, curr_event_idx);
		assert(band_kmer_to_offset(band_idx, curr_kmer_idx) == offset);

		uint8_t from = trace[band_idx][offset];
		if(from == FROM_D) {
		
			signalBuffer.push_back(r.events[curr_event_idx].mean);
			
			//if this query position is a match on the reference, use the kmer rank on the reference
			//accounts for basecalling inaccuracies under high analogue concentration
			if (signalBuffer.size() > 0 and r.queryToRef.count(curr_kmer_idx) > 0){
			
				unsigned int posOnRef = r.queryToRef.at(curr_kmer_idx);
				if (posOnRef < kmer_ranks_ref.size()){
					unsigned int kmer_rank_onRef = kmer_ranks_ref[posOnRef];
					cleanedRanks.push_back(kmer_rank_onRef);
					cleanedSignals.push_back(vectorMean(signalBuffer));
				}
			}
			signalBuffer.clear();
			
			curr_kmer_idx -= 1;
			curr_event_idx -= 1;
			curr_gap = 0;
			
		} else if(from == FROM_U) {
		
			signalBuffer.push_back(r.events[curr_event_idx].mean);

			curr_event_idx -= 1;
			curr_gap = 0;
		} else {
			curr_kmer_idx -= 1;
			curr_gap += 1;
			max_gap = std::max(curr_gap, max_gap);
		}   
	}
	std::reverse(r.eventAlignment.begin(), r.eventAlignment.end());

	//benchmarking
	//std::chrono::steady_clock::time_point tp4 = std::chrono::steady_clock::now();
	//std::cout << "backtrace: " << std::chrono::duration_cast<std::chrono::microseconds>(tp4 - tp3).count() << std::endl;

	// QC results
	double avg_log_emission = sum_emission / n_aligned_events;
	bool spanned = r.eventAlignment.front().second == 0 && r.eventAlignment.back().second == n_kmers - 1;
    
    	//Testing - print QCs
    	//std::cout << avg_log_emission << "\t" << spanned << "\t" << max_gap << "\t" << r.isReverse << std::endl;
    
	r.alignmentQCs.recordQCs(avg_log_emission, spanned, max_gap);
	
	//std::cerr << r.readID << " " << avg_log_emission << " " << max_gap << std::endl;

	std::pair<std::vector<double>, std::vector<unsigned int>> aligned_segmentation = std::make_pair(cleanedSignals, cleanedRanks);

	
	if(avg_log_emission < min_average_log_emission || !spanned || max_gap > max_gap_threshold ){
		r.eventAlignment.clear();
		return aligned_segmentation;
	}
	
	if ( cleanedSignals.size() < 1000 or cleanedRanks.size() < 1000){
		r.eventAlignment.clear();
		return aligned_segmentation;
	}
	
	return aligned_segmentation;

	//benchmarking
	//std::chrono::steady_clock::time_point tp5 = std::chrono::steady_clock::now();
	//std::cout << "calculate shift and scale: " << std::chrono::duration_cast<std::chrono::microseconds>(tp5 - tp4).count() << std::endl;
}


std::vector<double> quantileMedians(std::vector<double> &data, int nquantiles){

	auto endSlice = data.end();
	
	//uncomment to downsample to a fixed number of events
	/*
	unsigned int maxEvents = 100000;
	if (data.size() > maxEvents){
		endSlice = data.begin() + maxEvents;
	}
	*/
	
	std::vector<double> data_downsample(data.begin(), endSlice);
	std::sort(data_downsample.begin(), data_downsample.end());
	
	std::vector<double> quantileMedians;
	unsigned int n = data_downsample.size() / nquantiles;
	for (int i = 0; i < nquantiles; i++){
	
		double median  = data_downsample[ (i*n + (i+1)*n)/2 ];
		quantileMedians.push_back(median);
	}
	
	return quantileMedians;
}


std::pair<double, double> linear_regression(std::vector<double> x, std::vector<double> y){
//adapted from https://stackoverflow.com/questions/5083465/fast-efficient-least-squares-fit-algorithm-in-c

	assert(x.size() == y.size());
	
	double sum_x = 0., sum_x2 = 0., sum_y = 0., sum_xy = 0.;
	int n = y.size();
	
	for(int i = 0; i < n; i++){
		sum_x = sum_x + x[i];
		sum_x2 = sum_x2 + x[i]*x[i];
		sum_y = sum_y + y[i];
		sum_xy = sum_xy + x[i]*y[i];
	}
	
	//calculate coefficients
	double slope = (n * sum_xy - sum_x * sum_y)/(n * sum_x2 - sum_x * sum_x);
	double intercept = (sum_y - slope * sum_x)/n;
	
	//testing
	/*
	for(int i = 0; i < n; i++){
		std::cerr << x[i] << " " << y[i] << std::endl;
	}
	std::cerr << slope << " " << intercept << std::endl;
	std::cerr << "----------------------" << std::endl;	
	*/

	return std::make_pair(slope,intercept);
}


PoreParameters estimateScaling_quantiles(std::vector< double > &signal_means, std::string &sequence, std::vector<unsigned int> &kmer_ranks, bool useFitPoreModel ){

	PoreParameters s;

	std::vector<double> model_means;
	model_means.reserve(kmer_ranks.size());
	for ( unsigned int i = 0; i < kmer_ranks.size(); i ++ ){

		assert(kmer_ranks[i] < Pore_Substrate_Config.pore_model.size());
		
		std::pair<double,double> meanStd;
		if (useFitPoreModel){
			meanStd = Pore_Substrate_Config.unlabelled_model[kmer_ranks[i]];		
		}
		else{
			meanStd = Pore_Substrate_Config.pore_model[kmer_ranks[i]];
		}
		
		double kmer_mean = meanStd.first;
		model_means.push_back(kmer_mean);
	}

	std::vector<double> signal_quantiles = quantileMedians(signal_means, 10);
	std::vector<double> model_quantiles = quantileMedians(model_means, 10);

	std::pair<double, double> scalings = linear_regression(model_quantiles, signal_quantiles);
	
	s.shift = scalings.second;
	s.scale = scalings.first;
	
	return s;
}


void normaliseEvents( read &r, bool useFitPoreModel ){

	try{

		fast5_getSignal(r);
	}
	catch ( BadFast5Field &bf5 ){

		return;
	}
	
	event_table et = detect_events(&(r.raw)[0], (r.raw).size(), event_detection_defaults);
	assert(et.n > 0);
	
	r.events.reserve(et.n);
	unsigned int rawStart = 0;
	double mean = 0.;
	std::vector<double> event_means;
	for ( unsigned int i = 0; i < et.n; i++ ){

		if (et.event[i].mean > 0.) {

			if (i > 0){

				//build the previous event
				event e;
				e.mean = mean;
				event_means.push_back(mean);
				for (unsigned int j = rawStart; j <= std::min(et.event[i].start-1, r.raw.size()-1); j++){
				
					e.raw.push_back(r.raw[j]);
				}
				r.events.push_back(e);				
				
				//save stats for the next event
				mean = et.event[i].mean;
				rawStart = et.event[i].start;
			}
		}
	}
	free(et.event);
	
	// Precompute k-mer ranks for rescaling and banded alignment - query sequence
	size_t k = Pore_Substrate_Config.kmer_len;
	size_t n_kmers = r.basecall.size() - k + 1;
	std::vector<unsigned int> kmer_ranks_query(n_kmers);
	for(size_t i = 0; i < n_kmers; i++) {
		std::string kmer = r.basecall.substr(i, k);
		kmer_ranks_query[i] = kmer2index(kmer, k);
	}
	
	// Precompute k-mer ranks for rescaling and banded alignment - reference sequence
	n_kmers = r.referenceSeqMappedTo.size() - k + 1;
	std::vector<unsigned int> kmer_ranks_ref(n_kmers);
	for(size_t i = 0; i < n_kmers; i++) {
		std::string kmer = r.referenceSeqMappedTo.substr(i, k);
		kmer_ranks_ref[i] = kmer2index(kmer, k);
	}

	//normalise by quantile scaling by comparing the raw signal against the reference sequence
	r.scalings = estimateScaling_quantiles( event_means, r.referenceSeqMappedTo, kmer_ranks_ref, useFitPoreModel );

	// Rough alignment of signals to query sequence
	std::pair<std::vector<double>, std::vector<unsigned int>> segmentation = adaptive_banded_simple_event_align(r, kmer_ranks_query, kmer_ranks_ref, useFitPoreModel);

	//fine tune scaling parameters
	r.scalings = estimateScaling_theilSen(segmentation.first, segmentation.second, r.scalings, useFitPoreModel );
	
	//fail the read if it fails scaling refminement
	if (r.scalings.shift == -1.) r.eventAlignment.clear();

	r.scalings.eventsPerBase = (double) et.n / (double) (r.basecall.size() - k);
}
