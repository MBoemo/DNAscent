//----------------------------------------------------------
// Copyright 2017 University of Oxford
// Written by Michael A. Boemo (michael.boemo@path.ox.ac.uk)
// This software is licensed under GPL-3.0.  You should have
// received a copy of the license with this software.  If
// not, please Email the author.
//----------------------------------------------------------

#include <iterator>
#include "error_handling.h"



#define _USE_MATH_DEFINES


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

	for ( int i = 0; i < A.size(); i++ ){

		result.push_back( A[i].back() );
	}
	return result;
}


double fisherRaoMetric( double mu1, double stdv1, double mu2, double stdv2 ){

	return ( pow( mu1 - mu2, 2.0 ) + 2*pow( stdv1 - stdv2, 2.0 ) )/pow( stdv2, 2.0 );

}


std::vector< std::pair< int, std::string > > matchWarping( std::vector< double > &raw, 
							      std::vector< double > &raw_stdv,
						  	      std::string &basecall ){

	unsigned int numOfRaw = raw.size();
	unsigned int numOf5mers = basecall.size() - 4;

	std::vector< std::pair< int, std::string > > event5merPairs;

	/*allocate the dynamic time warping lattice */
	std::vector< std::vector< double > > dtw( numOfRaw, std::vector< double >( numOf5mers, std::numeric_limits< double >::max() ) );

	/*INITIALISATION */
	dtw[0][0] = 0.0;
	double mu = FiveMer_model[basecall.substr(1,5)].first;
	double stdv = FiveMer_model[basecall.substr(1,5)].second;
	dtw[1][1] = fisherRaoMetric( mu, stdv, raw[1], raw_stdv[1] );

	/*RECURSION: fill in the dynamic time warping lattice */
	for ( int row = 1; row < numOfRaw; row++ ){

		for ( int col = 2; col < numOf5mers; col++ ){

			mu = FiveMer_model[basecall.substr(col, 5)].first;
			stdv = FiveMer_model[basecall.substr(col, 5)].second;
			dtw[row][col] =  fisherRaoMetric( mu, stdv, raw[row], raw_stdv[row] ) + std::min( dtw[row - 1][col], std::min(dtw[row - 1][col - 1], dtw[row - 1][col - 2] ) );	
		}
	}

	/*TERMINATION: calculate the optimal warping path */
	int col = numOf5mers - 1;
	int row = numOfRaw - 1;

	while ( row > 0 and col > 1 ){

		event5merPairs.push_back( std::make_pair( row, basecall.substr(col, 5) ) );

		std::vector< double > mCand = { dtw[row - 1][col], dtw[row - 1][col - 1], dtw[row - 1][col - 2] };

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

	}
	std::reverse( event5merPairs.begin(), event5merPairs.end() );

	return event5merPairs;
}


std::vector< double > roughRescale( std::vector< double > means, std::string &basecalls ){

	unsigned int numOfFiveMers = basecalls.size() - 4;

	/*get a rough estimate for shift */
	double event_sum = 0.0;
	for ( unsigned int i = 0; i < means.size(); i++ ){

		event_sum += means[i];
	}

	double fiveMer_sum = 0.0;
	double fiveMer_sum_sq = 0.0;
	for ( unsigned int i = 0; i < numOfFiveMers; i ++ ){

		double fiveMer_mean = FiveMer_model[basecalls.substr(i, 5)].first;
		fiveMer_sum += fiveMer_mean;
		fiveMer_sum_sq = pow( fiveMer_mean, 2.0 );
	}

	double shift = event_sum / means.size() - fiveMer_sum / numOfFiveMers;

	/*rescale using shift */
	for ( unsigned int i = 0; i < means.size(); i++ ){

		means[i] = means[i] - shift;
	}
	return means;
}


std::vector< double > normaliseEvents( read r ){

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
	std::vector< double > rough_mu = roughRescale( events_mu, r.basecalls );

	/*align 5mers to events using the basecall */
	std::vector< std::pair< int, std::string > > event5merPairs = matchWarping( rough_mu, events_stdv, r.basecalls );

	/*calculate shift and scale */
	std::vector< std::vector< double > > A( 2, std::vector< double >( 2, 0.0 ) );
	std::vector< double > b( 2, 0.0 );

	for ( auto event = event5merPairs.begin(); event < event5merPairs.end(); event++ ){

		double event_mean = events_mu[ event -> first ];
		std::string fiveMer = event -> second;
		double model_mean = FiveMer_model[fiveMer].first;
		double model_stdv = FiveMer_model[fiveMer].second;

		A[0][0] += 1.0/pow( model_stdv, 2 );
		A[0][1] += 1.0/pow( model_stdv, 2 ) * model_mean;
		A[1][1] += 1.0/pow( model_stdv, 2 ) * pow( model_mean, 2 );

		b[0] += 1.0/pow( model_stdv, 2 ) * event_mean;
		b[1] += 1.0/pow( model_stdv, 2 ) * event_mean * model_mean;
	}

	/*use the symmetry of A */
	A[1][0] = A[0][1];

	/*solve the linear system and pull out values for shift and scale */
	std::vector< double > solution = solveLinearSystem( A, b );

	double shift = solution[0];
	double scale = solution[1];

	/*normalise event means for shift and scale */
	std::vector< double > normalisedEvents;
	normalisedEvents.reserve( event5merPairs.size() );

	double normalisedEventMean;
	for ( unsigned int i = 0 ; i < event5merPairs.size(); i++ ){

		normalisedEventMean = ( events_mu[i] - shift )/scale;

		if ( normalisedEventMean > 50.0 and normalisedEventMean < 130 ){

			normalisedEvents.push_back( normalisedEventMean );
		}
	}
	return normalisedEvents;
}


std::map< std::string, std::vector< std::vector< double > > > segmentEvents( std::string filename_foh, int threads ){

	std::map< std::string, std::vector< std::vector< double > > > scaled_trainingData;

	/*import the training data from the .foh file */
	std::vector< std::pair< std::string, std::vector< read > > > trainingData = import_foh( filename_foh );

	/*for each .foh entry, normalise the reads and add them to a map keyed by the .foh name */
	std::cout << "Normalising for shift and scale..." << std::endl;
	int prog = 0;
	int total = trainingData.size();

	#pragma omp parallel for default(none) shared(prog, total, scaled_trainingData, trainingData) num_threads(threads)
	for ( auto it = trainingData.begin(); it < trainingData.end(); it++ ){

		std::vector< std::vector< double > > parallelHolder;

		/*for each read grouped under this foh name */
		for ( auto r = (it -> second).begin(); r < (it -> second).end(); r++ ){

			std::vector< double > normalisedEvents = normaliseEvents( *r );

			parallelHolder.push_back( normalisedEvents );
		}

		#pragma omp atomic
		prog++;

		#pragma omp critical
		{	
			scaled_trainingData[ it -> first] = parallelHolder;
			displayProgress( prog, total );
		}
	}
	displayProgress( total, total );
	std::cout << std::endl << "Done." << std::endl;
	return scaled_trainingData;
}
