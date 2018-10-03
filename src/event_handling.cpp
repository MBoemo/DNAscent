//----------------------------------------------------------
// Copyright 2017 University of Oxford
// Written by Michael A. Boemo (michael.boemo@path.ox.ac.uk)
// This software is licensed under GPL-3.0.  You should have
// received a copy of the license with this software.  If
// not, please Email the author.
//----------------------------------------------------------

#include <iterator>
#include <algorithm>
#include <math.h>
#include "../Penthus/src/probability.h"
#include "error_handling.h"
#include "event_handling.h"

#define _USE_MATH_DEFINES

extern std::map< std::string, std::pair< double, double > > BrdU_model_full, SixMer_model;

std::vector< double > solveLinearSystem( std::vector< std::vector< double > > A, std::vector< double > b ){
/*crude but functional algorithm that solves the linear system A*x = b by building an augmented matrix and transforming to reduced row echelon form */
	
	if ( A.size() != b.size() ) throw MismatchedDimensions();

	int lead = 0;

	for ( unsigned int i = 0; i < A.size(); i++ ){

		A[i].push_back( b[i] );
	}

	int nRows = A.size() - 1;
	int nCols = A[0].size() - 1;

	for ( int r = 0; r <= nRows; ++r ){

		if ( lead > nCols ){
			break;
		}

		int i = r;

		while ( A[i][lead] == 0 ){
			
			++i;

			if ( i > nRows ){
			
				i = r;
				++lead;
				if ( lead > nCols ){
					goto stop;
				}
			}
		}

		std::vector< std::vector< double > > B = A;
		A[i] = B[r];
		A[r] = B[i];

		double normaliser = A[r][lead];
		if ( A[r][lead] != 0 ){
			
			for ( int c = 0; c <= nCols; ++c ){

				A[r][c] /= normaliser;
			}
		}

		for ( int i = 0; i <= nRows; ++i ){

			if ( i != r ){
				normaliser = A[i][lead];
				for ( int c = 0; c <= nCols; ++c ){ 
					A[i][c] -= normaliser*A[r][c];
				}
			}
		}
		lead++;
	}

	stop:
	std::vector< double > result;
	result.reserve( A.size() );

	for ( unsigned int i = 0; i < A.size(); i++ ){

		result.push_back( A[i].back() );
	}
	return result;
}


double fisherRaoMetric( double mu1, double stdv1, double mu2, double stdv2 ){
/*computes the length of the geodesic between N(mu1,stdv1) and N(mu2,stdv2) using the Fisher-Rao metric as a distance */

	double F = sqrt( ( pow( mu1 - mu2, 2.0 ) + 2*pow( stdv1 - stdv2, 2.0 ) )*( pow( mu1 - mu2, 2.0 ) + 2*pow( stdv1 + stdv2, 2.0 ) ) );

	return sqrt(2)*log( (F + pow( mu1 - mu2, 2.0 ) + 2*( pow( stdv1, 2.0 ) + pow( stdv2, 2.0 ) ) ) / ( 4 * stdv1 * stdv2 ) );
}


double manhattanMetric( double mu1, double stdv1, double mu2, double stdv2 ){
/*computes the length of the geodesic between N(mu1,stdv1) and N(mu2,stdv2) using the Fisher-Rao metric as a distance */

	return fabs(mu1-mu2);
}


std::vector< std::pair< unsigned int, unsigned int > > matchWarping( std::vector< double > &raw, std::vector< double > &raw_stdv, std::string &basecall ){
/*use dynamic time warping to calculate an alignment between the raw signal and the basecall */

	unsigned int numOfRaw = raw.size();
	unsigned int numOf5mers = basecall.size() - 5;

	std::vector< std::pair< unsigned int, unsigned int > > eventSeqLocPairs;

	/*allocate the dynamic time warping lattice */
	std::vector< std::vector< double > > dtw( numOfRaw, std::vector< double >( numOf5mers, std::numeric_limits< double >::max() ) );

	/*INITIALISATION */
	dtw[0][0] = 0.0;
	dtw[1][1] = 0.0;
	double mu = SixMer_model[basecall.substr(1,6)].first;
	double stdv = SixMer_model[basecall.substr(1,6)].second;
	dtw[1][1] = manhattanMetric( mu, stdv, raw[1], raw_stdv[1] );

	/*RECURSION: fill in the dynamic time warping lattice */
	for ( unsigned int row = 1; row < numOfRaw; row++ ){

		for ( unsigned int col = 2; col < numOf5mers; col++ ){

			mu = SixMer_model[basecall.substr(col, 6)].first;
			stdv = SixMer_model[basecall.substr(col, 6)].second;
			dtw[row][col] =  manhattanMetric( mu, stdv, raw[row], raw_stdv[row] ) + std::min( dtw[row - 1][col], std::min(dtw[row - 1][col - 1], dtw[row - 1][col - 2] ) );	
		}
	}

	/*TRACEBACK: calculate the optimal warping path */
	int col = numOf5mers - 1;
	int row = numOfRaw - 1;
	eventSeqLocPairs.push_back( std::make_pair( row, col ) );

	while ( row > 0 and col > 0 ){

		std::vector< double > mCand;

		if (col == 1){
			mCand = { dtw[row - 1][col], dtw[row - 1][col - 1] };
		}
		else {
			mCand = { dtw[row - 1][col], dtw[row - 1][col - 1], dtw[row - 1][col - 2] };
		}

		int m = std::min_element( mCand.begin(), mCand.end() ) - mCand.begin();

		if ( m == 0 ){
			row--;
		}
		else if ( m == 1 ){
			row--;
			col--;
		}
		else if ( m == 2 ){
			row--;
			col-=2;
		}
		else{
			std::cout << "Exiting with error.  Out of bounds error in dynmaic time warping." << std::endl;
			exit(EXIT_FAILURE);
		}

		eventSeqLocPairs.push_back( std::make_pair( row, col ) );
	}
	std::reverse( eventSeqLocPairs.begin(), eventSeqLocPairs.end() );

	return eventSeqLocPairs;
}


std::vector< std::pair< unsigned int, unsigned int > > matchWarping_band( std::vector< double > &raw, std::vector< double > &raw_stdv, std::string &basecall ){
/*use dynamic time warping to calculate an alignment between the raw signal and the basecall */

	int numOfRaw = raw.size();
	int numOf5mers = basecall.size() - 5;
	int bandWidth = 500;

	std::vector< std::pair< unsigned int, unsigned int > > eventSeqLocPairs;

	/*allocate the dynamic time warping lattice */
	std::vector< std::vector< double > > dtw( numOfRaw, std::vector< double >( numOf5mers, std::numeric_limits< double >::max() ) );

	/*INITIALISATION */
	dtw[0][0] = 0.0;
	dtw[1][1] = 0.0;
	double mu = SixMer_model[basecall.substr(1,6)].first;
	double stdv = SixMer_model[basecall.substr(1,6)].second;
	dtw[1][1] = manhattanMetric( mu, stdv, raw[1], raw_stdv[1] );

	/*RECURSION: fill in the dynamic time warping lattice */
	for ( unsigned int col = 2; col < numOf5mers; col++ ){//basecall

		int bandLower = col*numOfRaw/numOf5mers - bandWidth;
		int bandUpper = col*numOfRaw/numOf5mers + bandWidth;

		bandLower = std::max(1, bandLower);
		bandUpper = std::min(numOfRaw, bandUpper);
		
		for ( unsigned int row = bandLower; row < bandUpper; row++ ){//raw

			mu = SixMer_model[basecall.substr(col, 6)].first;
			stdv = SixMer_model[basecall.substr(col, 6)].second;
			dtw[row][col] =  manhattanMetric( mu, stdv, raw[row], raw_stdv[row] ) + std::min( dtw[row - 1][col], std::min(dtw[row - 1][col - 1], dtw[row - 1][col - 2] ) );	
		}
	}

	/*TRACEBACK: calculate the optimal warping path */
	int col = numOf5mers - 1;
	int row = numOfRaw - 1;
	eventSeqLocPairs.push_back( std::make_pair( row, col ) );

	while ( row > 0 and col > 0 ){

		std::vector< double > mCand;

		if (col == 1){
			mCand = { dtw[row - 1][col], dtw[row - 1][col - 1] };
		}
		else {
			mCand = { dtw[row - 1][col], dtw[row - 1][col - 1], dtw[row - 1][col - 2] };
		}

		int m = std::min_element( mCand.begin(), mCand.end() ) - mCand.begin();

		if ( m == 0 ){
			row--;
		}
		else if ( m == 1 ){
			row--;
			col--;
		}
		else if ( m == 2 ){
			row--;
			col-=2;
		}
		else{
			std::cout << "Exiting with error.  Out of bounds error in dynmaic time warping." << std::endl;
			exit(EXIT_FAILURE);
		}

		eventSeqLocPairs.push_back( std::make_pair( row, col ) );
	}
	std::reverse( eventSeqLocPairs.begin(), eventSeqLocPairs.end() );

	return eventSeqLocPairs;
}


//start: adapted from nanopolish

inline float logProbabilityMatch(std::string sixMer, double x){

	if (sixMer.find('T') != std::string::npos and BrdU_model_full.count(sixMer) > 0){

		double mu = SixMer_model.at(sixMer).first;
		double sigma = SixMer_model.at(sixMer).second;
		double thymProb = ( 1.0/sqrt( 2.0*pow( sigma, 2.0 )*M_PI ) )*exp( -pow( x - mu , 2.0 )/( 2.0*pow( sigma, 2.0 ) ) );

		mu = BrdU_model_full.at(sixMer).first;
		sigma = BrdU_model_full.at(sixMer).second;
		double brduProb = ( 1.0/sqrt( 2.0*pow( sigma, 2.0 )*M_PI ) )*exp( -pow( x - mu , 2.0 )/( 2.0*pow( sigma, 2.0 ) ) );

		if ( thymProb > brduProb ) return eln(thymProb);
		else return eln(brduProb);
	}
	else{
		double mu = SixMer_model.at(sixMer).first;
		double sigma = SixMer_model.at(sixMer).second;

		double prob = ( 1.0/sqrt( 2.0*pow( sigma, 2.0 )*M_PI ) )*exp( -pow( x - mu , 2.0 )/( 2.0*pow( sigma, 2.0 ) ) );

		return eln(prob);
	}
}

#define event_kmer_to_band(ei, ki) (ei + 1) + (ki + 1)
#define band_event_to_offset(bi, ei) band_lower_left[bi].event_idx - (ei)
#define band_kmer_to_offset(bi, ki) (ki) - band_lower_left[bi].kmer_idx
#define is_offset_valid(offset) (offset) >= 0 && (offset) < bandwidth
#define event_at_offset(bi, offset) band_lower_left[(bi)].event_idx - (offset)
#define kmer_at_offset(bi, offset) band_lower_left[(bi)].kmer_idx + (offset)
#define move_down(curr_band) { curr_band.event_idx + 1, curr_band.kmer_idx }
#define move_right(curr_band) { curr_band.event_idx, curr_band.kmer_idx + 1 }

std::vector< std::pair< unsigned int, unsigned int > > adaptive_banded_simple_event_align( std::vector< double > &raw, const std::string &sequence ){

	size_t strand_idx = 0;
	size_t k = 6;
	//const Alphabet* alphabet = pore_model.pmalphabet;
	size_t n_events = raw.size();
	size_t n_kmers = sequence.size() - k + 1;

	// backtrack markers
	const uint8_t FROM_D = 0;
	const uint8_t FROM_U = 1;
	const uint8_t FROM_L = 2;
 
	// qc
	double min_average_log_emission = -10;//-5.0;
	int max_gap_threshold = 50;

	// banding
	//int bandwidth = 100;
	//scale bandwidth with sequence length
	int bandwidth;
	if ( sequence.length() < 10000 ){
	
		bandwidth = 100;
	}
	else {

		bandwidth = (int) sequence.length() / 100;
		if (bandwidth % 2 == 1) bandwidth++;//make it even
	}

	int half_bandwidth = bandwidth / 2;
 
	// transition penalties
	double events_per_kmer = (double)n_events / n_kmers;
	double p_stay = 1 - (1 / (events_per_kmer + 1));

	// setting a tiny skip penalty helps keep the true alignment within the adaptive band
	// this was empirically determined
	double epsilon = 1e-10;
	double lp_skip = log(epsilon);
	double lp_stay = log(p_stay);
	double lp_step = log(1.0 - exp(lp_skip) - exp(lp_stay));
	double lp_trim = log(0.01);
 
	// dp matrix
	size_t n_rows = n_events + 1;
	size_t n_cols = n_kmers + 1;
	size_t n_bands = n_rows + n_cols;
 
	// Initialize

	// Precompute k-mer ranks to avoid doing this in the inner loop
	//std::vector<size_t> kmer_ranks(n_kmers);
	//for(size_t i = 0; i < n_kmers; ++i) {
	//	kmer_ranks[i] = sixMerRank_nanopolish(sequence.substr(i, k).c_str());
	//}

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
	for(int band_idx = 2; band_idx < n_bands; ++band_idx) {
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
			int event_idx = event_at_offset(band_idx, trim_offset);
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

			//size_t kmer_rank = kmer_ranks[kmer_idx];
			std::string sixMer = sequence.substr(kmer_idx, k);
 
			int offset_up   = band_event_to_offset(band_idx - 1, event_idx - 1); 
			int offset_left = band_kmer_to_offset(band_idx - 1, kmer_idx - 1);
			int offset_diag = band_kmer_to_offset(band_idx - 2, kmer_idx - 1);

			float up   = is_offset_valid(offset_up)   ? bands[band_idx - 1][offset_up]   : -INFINITY;
			float left = is_offset_valid(offset_left) ? bands[band_idx - 1][offset_left] : -INFINITY;
			float diag = is_offset_valid(offset_diag) ? bands[band_idx - 2][offset_diag] : -INFINITY;
 
			//float lp_emission = log_probability_match_r9(read, pore_model, kmer_rank, event_idx, strand_idx);
			float lp_emission = logProbabilityMatch(sixMer, raw[event_idx]);

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

	//
	// Backtrack to compute alignment
	//
	double sum_emission = 0;
	double eventDiffs = 0.0;
	double n_aligned_events = 0;
	//std::vector<AlignedPair> out;
	std::vector< std::pair< unsigned int, unsigned int > > eventSeqLocPairs;
    
	float max_score = -INFINITY;
	int curr_event_idx = 0;
	int curr_kmer_idx = n_kmers -1;

	// Find best score between an event and the last k-mer. after trimming the remaining evnets
	for(int event_idx = 0; event_idx < n_events; ++event_idx) {
		int band_idx = event_kmer_to_band(event_idx, curr_kmer_idx);
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

	int curr_gap = 0;
	int max_gap = 0;
	while(curr_kmer_idx >= 0 && curr_event_idx >= 0) {
        
		// emit alignment
		eventSeqLocPairs.push_back(std::make_pair(curr_event_idx, curr_kmer_idx));

		// qc stats
		//size_t kmer_rank = sixMerRank_nanopolish(sequence.substr(curr_kmer_idx, k).c_str());
		//sum_emission += log_probability_match_r9(read, pore_model, kmer_rank, curr_event_idx, strand_idx);
		std::string sixMer = sequence.substr(curr_kmer_idx, k);
		sum_emission += logProbabilityMatch(sixMer, raw[curr_event_idx]);
		eventDiffs += SixMer_model.at(sixMer).first - raw[curr_event_idx];

		n_aligned_events += 1;

		int band_idx = event_kmer_to_band(curr_event_idx, curr_kmer_idx);
		int offset = band_event_to_offset(band_idx, curr_event_idx);
		assert(band_kmer_to_offset(band_idx, curr_kmer_idx) == offset);

		uint8_t from = trace[band_idx][offset];  //PROBLEM HERE
		if(from == FROM_D) {
			curr_kmer_idx -= 1;
			curr_event_idx -= 1;
			curr_gap = 0;
		} else if(from == FROM_U) {
			curr_event_idx -= 1;
			curr_gap = 0;
		} else {
			curr_kmer_idx -= 1;
			curr_gap += 1;
			max_gap = std::max(curr_gap, max_gap);
		}   
	}
	std::reverse(eventSeqLocPairs.begin(), eventSeqLocPairs.end());
    
	// QC results
	double avg_log_emission = sum_emission / n_aligned_events;
	double avg_eventDiffs = eventDiffs / n_aligned_events;
	bool spanned = eventSeqLocPairs.front().second == 0 && eventSeqLocPairs.back().second == n_kmers - 1;
    
	bool failed = false;
	if(avg_log_emission < min_average_log_emission || !spanned || max_gap > max_gap_threshold) {
		
		failed = true;		
		eventSeqLocPairs.clear();
    		//fprintf(stderr, "ada\t\t%s\t%s\t%.2lf\t%zu\t%.2lf\t%.2lf\t%d\t%d\t%d\n", failed ? "FAILED" : "OK",spanned ? "SPAN" : "NOTS", events_per_kmer, sequence.size(), avg_log_emission, avg_eventDiffs, curr_event_idx, max_gap, fills);
	}


	return eventSeqLocPairs;
}
//end:from nanopolish


std::vector< double > roughRescale( std::vector< double > means, std::string &basecall ){

	unsigned int numOfSixMers = basecall.size() - 5;

	/*get a rough estimate for shift */
	double event_sum = 0.0;
	for ( unsigned int i = 0; i < means.size(); i++ ){

		event_sum += means[i];
	}

	double sixMer_sum = 0.0;
	double sixMer_sq_sum = 0.0;
	for ( unsigned int i = 0; i < numOfSixMers; i ++ ){

		double sixMer_mean = SixMer_model[basecall.substr(i, 6)].first;
		sixMer_sum += sixMer_mean;
		sixMer_sq_sum += pow( sixMer_mean, 2.0 );
	}

	double shift = event_sum / means.size() - sixMer_sum / numOfSixMers;

	/*get a rough estimate for scale */
	double event_sq_sum = 0.0;
	for ( unsigned int i = 0; i < means.size(); i++ ){

		event_sq_sum += pow( means[i] - shift, 2.0 );
	}

	double scale = (event_sq_sum / means.size()) / (sixMer_sq_sum / numOfSixMers ); 

	/*reestimate event means adjusting for shift and scale */
	for ( unsigned int i = 0; i < means.size(); i++ ){

		means[i] = (means[i] - shift) / scale;
	}
	return means;
}


void normaliseEvents( read &r ){

	/*allocate some space for event detection */
	event_s *c_events = (event_s *)calloc( (r.raw).size(), sizeof( event_s) );

	/*we can trim and segment the raw signal if needed - leave these commented out for now */
	//unsigned int rawStart, rawEnd;
	//trim_and_segment_raw( &(r.raw)[0], (r.raw).size(), &rawStart, &rawEnd );

	/*use Scrappie to segment events based on the raw signal */
	size_t numOfDetectedEvents;
	detect_events( &(r.raw)[0], (r.raw).size(), event_detection_defaults, c_events, &numOfDetectedEvents );

	/*we only care about the mean and stdv, so pull these out so they're easier to work with */
	std::vector< double > events_mu;
	events_mu.reserve( numOfDetectedEvents );
	std::vector< double > events_stdv;
	events_stdv.reserve( numOfDetectedEvents );

	for ( int i = 0; c_events[i].mean != 0; i++ ){

		events_mu.push_back( c_events[i].mean );
		events_stdv.push_back( c_events[i].stdv );
	}
	free(c_events);

	/*rough calculation of shift and scale so that we can align events */
	std::vector< double > rough_mu = roughRescale( events_mu, r.basecall );

	r.normalisedEvents = rough_mu;

	/*align 5mers to events using the basecall */
	r.eventAlignment = adaptive_banded_simple_event_align(rough_mu, r.basecall);

	double normalisedEventMean;
	double alignmentScore = 0.0;
	int positionAlignedTo;
	std::string sixMerAlignedTo;
	int numEventsAdded = 0;
	for ( unsigned int i = 0 ; i < (r.eventAlignment).size(); i++ ){

		normalisedEventMean = rough_mu[(r.eventAlignment)[i].first];

		positionAlignedTo = (r.eventAlignment)[i].second;
		sixMerAlignedTo = (r.basecall).substr(positionAlignedTo, 6);

		alignmentScore += normalisedEventMean - SixMer_model[sixMerAlignedTo].first;
		numEventsAdded++;

		//std::cout << positionAlignedTo << '\t' << normalisedEventMean << '\t' << sixMerAlignedTo << '\t' << SixMer_model[sixMerAlignedTo].first << '\t' << SixMer_model[sixMerAlignedTo].second << std::endl;
	}
	//std::cout << "-------------" << std::endl;
	r.qualityScore = alignmentScore / (double) numEventsAdded;
}
