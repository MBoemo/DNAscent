//----------------------------------------------------------
// Copyright 2017 University of Oxford
// Written by Michael A. Boemo (michael.boemo@path.ox.ac.uk)
//----------------------------------------------------------


#include "common.h"


int argMin( std::vector< double > &vec ){

	double smallest = vec[0];
	int index;

	for ( unsigned int i = 1; i < vec.size(); i++ ){

		if ( vec[i] < smallest ){

			smallest = vec[i];
			index = i;

		}
	
	}

	return index;

}


std::vector< double > generateSignal( std::string &reference, std::map< std::string, std::pair< double, double > > &baseModel ){

	std::vector< double > generatedSignal;
	generatedSignal.reserve( reference.length() - 6 );

	for ( unsigned int i = 0; i < reference.length() - 6; i++ ){

		generatedSignal.push_back( ( baseModel.at( reference.substr( i, 6 ) ) ).first );

	}

	return generatedSignal;

}


std::vector< std::vector< double > > filterEvents( std::string &reference, std::map< std::string, std::pair< double, double > > &baseModel, std::vector< std::vector< double > > &events ){

	std::vector< double > generatedSignal = generateSignal( reference, baseModel );

	std::vector< std::vector< double > > filteredEvents;

	for ( unsigned int i = 0; i < events.size(); i++ ){

		std::pair< int, int > subEventBounds = subsequenceDynamicTimewarping( generatedSignal, events[ i ] );

		if (subEventBounds.second != 0 ){

			std::vector< double > subEvent( events[ i ].begin() + subEventBounds.first, events[ i ].begin() + subEventBounds.second );

			filteredEvents.push_back( subEvent );

		}

	}

	return filteredEvents;

}


std::pair< int, int > subsequenceDynamicTimewarping( std::vector< double > &x, std::vector< double > &y ){

	/*allocate the dynamic time warping lattice */
	std::vector< std::vector< double > > dtw( x.size(), std::vector< double >( y.size(), 0.0 ) );

	/*initialisation: set the first row and first column of the dynamic time warping lattice */
	double runningSum = 0.0;
	for ( int i = 0; i < x.size(); i++ ){

		runningSum += std::abs( x[ i ] - y[ 0 ] );
		dtw[ i ][ 0 ] = runningSum;

	}

	for ( int i = 0; i < y.size(); i++ ){

		dtw[ 0 ][ i ] = std::abs( x[ 0 ] - y[ i ] );

	}

	/*recursion: fill in the dynamic time warping lattice */
	for ( int i = 1; i < x.size(); i++ ){
		for ( int j = 1; j < y.size(); j++ ){

			dtw[ i ][ j ] = std::abs( x[ i ] - y[ j ] ) + std::min( dtw[ i - 1 ][ j ], std::min( dtw[ i ][ j - 1 ], dtw[ i - 1 ][ j - 1 ] ) );

		}
	
	}

	/*termination: calculate the optimal warping path */
	int j = argMin( dtw.back() );
	int i = x.size() - 1;
	std::vector< std::pair< int, int > > path;
	while ( i != 0 ){

		path.push_back( std::make_pair( i, j ) );

		std::vector< double > mCand = { dtw[ i - 1 ][ j - 1 ], dtw[ i ][ j - 1 ], dtw[ i - 1 ][ j ] };

		int m = std::min_element( mCand.begin(), mCand.end() ) - mCand.begin();

		if ( m == 0 ){
			j--;
			i--;
		}
		else if ( m == 1 ){
			i--;
		}
		else if ( m == 2 ){
			j--;
		}
		else{
			std::cout << "Exiting with error.  Out of bounds error in dynmaic time warping." << std::endl;
			exit(EXIT_FAILURE);
		}

		if ( j == 0 and i > 0 ){

			return std::make_pair( 0, 0 );

		}

	}

	return std::make_pair( ( path.back() ).second, ( path.front() ).second );

}
