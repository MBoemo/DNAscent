//----------------------------------------------------------
// Copyright 2017 University of Oxford
// Written by Michael A. Boemo (michael.boemo@path.ox.ac.uk)
// This software is licensed under GPL-3.0.  You should have
// received a copy of the license with this software.  If
// not, please Email the author.
//----------------------------------------------------------

#include "scrappie/event_detection.h"
#include "scrappie/scrappie_common.h"
#include "poreModels.h"
#include "common.h"
#include "data_IO.h"
#include <math.h>
#include <stdlib.h>
#include <assert.h>


#define _USE_MATH_DEFINES


std::vector< double > reducedRowEchelonForm( std::vector< std::vector< double > > A, std::vector< double >b ){
/*crude but functional algorithm that solves the linear system A*x = b */
	
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


std::vector< std::pair< double, std::string > > matchWarping( std::vector< double > &raw, 
							      std::vector< double > &raw_stdv,
						  	      std::string &basecall ){

	unsigned int numOfRaw = raw.size();
	unsigned int numOf5mers = basecall.size() - 4;

	std::vector< std::pair< double, std::string > > event5merPairs;

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
			dtw[row][col] =  fisherRaoMetric( mu, stdv, raw[row], raw_stdv[row] ) + std::min( dtw[row - 1][col], std::min(dtw[row - 1][col - 1], 1.5*dtw[row - 1][col - 2] ) );	
		}
	}

	/*TERMINATION: calculate the optimal warping path */
	int col = numOf5mers - 1;
	int row = numOfRaw - 1;

	while ( row > 0 and col > 1 ){

		event5merPairs.push_back( std::make_pair( raw[row], basecall.substr(col, 5) ) );

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

	/*align 5mers to events using the basecall */
	std::vector< std::pair< double, std::string > > event5merPairs = matchWarping( events_mu, events_stdv, r.basecalls );


	/*normalise for shift/scale/drift */
	std::vector< std::vector< double > > A( 3, std::vector< double >( 3, 0.0 ) );
	std::vector< double > b( 3, 0.0 );

	for ( auto event = event5merPairs.begin(); event < event5merPairs.end(); event++ ){

		//build linear system

	}



	for ( unsigned int i = 0; i < event5merPairs.size(); i++ ){
	
		std::cout << event5merPairs[i].first << " " << event5merPairs[i].second << " " << FiveMer_model[event5merPairs[i].second].first << std::endl;
	}
	std::cout << "--------------------------------------------------------" << std::endl;
}


std::map< std::string, std::vector< std::vector< double > > > segmentEvents( std::string filename_foh ){

	std::map< std::string, std::vector< std::vector< double > > > scaled_trainingData;

	std::map< std::string, std::vector< read > > trainingData = import_foh( filename_foh );

	for ( auto it = trainingData.cbegin(); it != trainingData.cend(); it++ ){

		for ( auto r = (it -> second).begin(); r < (it -> second).end(); r++ ){

			scaled_trainingData[ it -> first].push_back( normaliseEvents( *r ) );
		}
	}
	return scaled_trainingData;
}
