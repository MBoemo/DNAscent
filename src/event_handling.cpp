//----------------------------------------------------------
// Copyright 2017 University of Oxford
// Written by Michael A. Boemo (michael.boemo@path.ox.ac.uk)
// This software is licensed under GPL-3.0.  You should have
// received a copy of the license with this software.  If
// not, please Email the author.
//----------------------------------------------------------

#include <iterator>
#include <math.h>
#include "error_handling.h"
#include "event_handling.h"

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
	unsigned int numOf5mers = basecall.size() - 4;

	std::vector< std::pair< unsigned int, unsigned int > > eventSeqLocPairs;

	/*allocate the dynamic time warping lattice */
	std::vector< std::vector< double > > dtw( numOfRaw, std::vector< double >( numOf5mers, std::numeric_limits< double >::max() ) );

	/*INITIALISATION */
	dtw[0][0] = 0.0;
	dtw[1][1] = 0.0;
	double mu = FiveMer_model[basecall.substr(1,5)].first;
	double stdv = FiveMer_model[basecall.substr(1,5)].second;
	dtw[1][1] = manhattanMetric( mu, stdv, raw[1], raw_stdv[1] );

	/*RECURSION: fill in the dynamic time warping lattice */
	for ( unsigned int row = 1; row < numOfRaw; row++ ){

		for ( unsigned int col = 2; col < numOf5mers; col++ ){

			mu = FiveMer_model[basecall.substr(col, 5)].first;
			stdv = FiveMer_model[basecall.substr(col, 5)].second;
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


std::vector< double > roughRescale( std::vector< double > means, std::string &basecalls ){

	unsigned int numOfFiveMers = basecalls.size() - 4;

	/*get a rough estimate for shift */
	double event_sum = 0.0;
	for ( unsigned int i = 0; i < means.size(); i++ ){

		event_sum += means[i];
	}

	double fiveMer_sum = 0.0;
	for ( unsigned int i = 0; i < numOfFiveMers; i ++ ){

		double fiveMer_mean = FiveMer_model[basecalls.substr(i, 5)].first;
		fiveMer_sum += fiveMer_mean;
	}

	double shift = event_sum / means.size() - fiveMer_sum / numOfFiveMers;

	/*rescale using shift */
	for ( unsigned int i = 0; i < means.size(); i++ ){

		means[i] = means[i] - shift;
	}
	return means;
}


eventDataForRead normaliseEvents( read &r ){

	eventDataForRead thisRead;

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
	thisRead.eventAlignment = matchWarping( rough_mu, events_stdv, r.basecalls );

	/*calculate shift and scale */
	std::vector< std::vector< double > > A( 2, std::vector< double >( 2, 0.0 ) );
	std::vector< double > b( 2, 0.0 );

	for ( auto event = (thisRead.eventAlignment).begin(); event < (thisRead.eventAlignment).end(); event++ ){

		double event_mean = events_mu[ event -> first ];
		std::string fiveMer = (r.basecalls).substr(event -> second, 5);
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
	(thisRead.normalisedEvents).reserve( (thisRead.eventAlignment).size() );

	double normalisedEventMean;
	double alignmentScore = 0.0;
	int positionAlignedTo;
	std::string fiveMerAlignedTo;
	int numEventsAdded = 0;
	for ( unsigned int i = 0 ; i < (thisRead.eventAlignment).size(); i++ ){

		normalisedEventMean = ( events_mu[i] - shift )/scale;

		positionAlignedTo = (thisRead.eventAlignment)[i].second;
		fiveMerAlignedTo = (r.basecalls).substr(positionAlignedTo,5);

		if ( positionAlignedTo >= (r.bounds_query).first and positionAlignedTo <= (r.bounds_query).second ){

			(thisRead.normalisedEvents).push_back( normalisedEventMean );
			alignmentScore += normalisedEventMean - FiveMer_model[fiveMerAlignedTo].first;
			numEventsAdded++;

			//std::cout << normalisedEventMean << '\t' << fiveMerAlignedTo << '\t' << FiveMer_model[fiveMerAlignedTo].first << '\t' << FiveMer_model[fiveMerAlignedTo].second << std::endl;

		}
	}
	thisRead.qualityScore = alignmentScore / (double) numEventsAdded;

	return thisRead;
}
