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


std::map< int, std::vector< int > > dynamicTimewarping( std::vector< double > &events, std::vector< double > &generatedSignal ){

	/*allocate the dynamic time warping lattice */
	std::vector< std::vector< double > > dtw( events.size(), std::vector< double >( generatedSignal.size(), 0.0 ) );

	/*initialisation */
	double runningSum = 0.0;
	for ( unsigned int i = 0; i < events.size(); i++ ){

		runningSum += std::abs( events[ i ] - generatedSignal[ 0 ] );
		dtw[ i ][ 0 ] = runningSum;

	}

	runningSum = 0.0;
	for ( unsigned int i = 0; i < generatedSignal.size(); i++ ){

		runningSum += std::abs( events[ 0 ] - generatedSignal[ i ] );
		dtw[ 0 ][ i ] = runningSum; 

	}

	/*recursion: fill in the dynamic time warping lattice */
	for ( int i = 1; i < events.size(); i++ ){
		for ( int j = 1; j < generatedSignal.size(); j++ ){

			if ( i == 1 or j == 1 ){
				dtw[ i ][ j ] = std::abs( events[ i ] - generatedSignal[ j ] ) + std::min( dtw[ i - 1 ][ j ], std::min( dtw[ i ][ j - 1 ], dtw[ i - 1 ][ j - 1 ] ) );			
			}
			else{
				dtw[ i ][ j ] = std::abs( events[ i ] - generatedSignal[ j ] ) + std::min( dtw[ i - 1 ][ j - 1 ], std::abs( events[ i - 1 ] - generatedSignal[ j ] ) + std::min( dtw[ i - 2 ][ j - 1 ], std::min( dtw[ i - 1 ][ j ], std::abs( events[ i - 1 ] - generatedSignal[ j ] ) + dtw[ i - 2 ][ j ] ) ) );

			}
		}
	
	}

	/*termination: calculate the optimal warping path */
	int j = generatedSignal.size() - 1;
	int i = events.size() - 1;
	std::map< int, std::vector< int > > readPosToEventPos;

	while ( i > 2 and j > 2 ){

		readPosToEventPos[ j ].push_back( i );
		std::cout << j << " " << i << std::endl;

		std::vector< double > mCand = { dtw[ i - 1 ][ j - 1 ], dtw[ i - 2 ][ j - 1 ], dtw[ i - 1 ][ j ], dtw[ i - 2 ][ j ] };

		int m = std::min_element( mCand.begin(), mCand.end() ) - mCand.begin();

		if ( m == 0 ){
			i--;
			j--;
		}
		else if ( m == 1 ){
			i-=2;
			j--;
		}
		else if ( m == 2 ){
			i--;
		}
		else if ( m == 3 ){
			i-=2;
		}
		else{
			std::cout << "Exiting with error.  Out of bounds error in dynmaic time warping." << std::endl;
			exit(EXIT_FAILURE);
		}

	}

	return readPosToEventPos;

}


std::pair< int, int > subsequenceDynamicTimewarping( std::vector< double > &shortSignal, std::vector< double > &longSignal ){

	/*allocate the dynamic time warping lattice */
	std::vector< std::vector< double > > dtw( shortSignal.size(), std::vector< double >( longSignal.size(), 0.0 ) );

	/*initialisation: set the first row and first column of the dynamic time warping lattice */
	double runningSum = 0.0;
	for ( int i = 0; i < shortSignal.size(); i++ ){

		runningSum += std::abs( shortSignal[ i ] - longSignal[ 0 ] );
		dtw[ i ][ 0 ] = runningSum;

	}

	for ( int i = 0; i < longSignal.size(); i++ ){

		dtw[ 0 ][ i ] = std::abs( shortSignal[ 0 ] - longSignal[ i ] );

	}

	/*recursion: fill in the dynamic time warping lattice */
	for ( int i = 1; i < shortSignal.size(); i++ ){
		for ( int j = 1; j < longSignal.size(); j++ ){

			dtw[ i ][ j ] = std::abs( shortSignal[ i ] - longSignal[ j ] ) + std::min( dtw[ i - 1 ][ j ], std::min( dtw[ i ][ j - 1 ], dtw[ i - 1 ][ j - 1 ] ) );

		}
	
	}

	/*termination: calculate the optimal warping path */
	int j = argMin( dtw.back() );
	int i = shortSignal.size() - 1;
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
			std::cout << "Exiting with error.  Out of bounds error in subsequence dynmaic time warping." << std::endl;
			exit(EXIT_FAILURE);
		}

		if ( j == 0 and i > 0 ){

			return std::make_pair( 0, 0 );

		}

	}

	return std::make_pair( ( path.back() ).second, ( path.front() ).second );

}
