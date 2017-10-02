//----------------------------------------------------------
// Copyright 2017 University of Oxford
// Written by Michael A. Boemo (michael.boemo@path.ox.ac.uk)
//----------------------------------------------------------


#include "common.h"


void displayProgress( int current, int total ){
/*poached from https://stackoverflow.com/questions/14539867/how-to-display-a-progress-indicator-in-pure-c-c-cout-printf */

	double progress = (double) current / (double) total;
	int barWidth = 70;

	if ( progress < 1.0 ){

		std::cout << "[";
		int pos = barWidth * progress;
		for (int i = 0; i < barWidth; ++i) {
			if (i < pos) std::cout << "=";
			else if (i == pos) std::cout << ">";
			else std::cout << " ";
		}
		std::cout << "] " << int(progress * 100.0) << " %\r";
		std::cout.flush();

	}

}


std::vector< std::string > split( std::string s, char delim ){

	std::stringstream ssString( s );
	std::vector< std::string > splitString;
	std::string entry;

	while ( std::getline( ssString, entry, delim ) ){

		splitString.push_back( entry );

	}

	return splitString;

}


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

		std::string sixMer = reference.substr( i, 6 );

		/*if the signal contains an N, use the median pA and let DTW find the indel.  Otherwise, use the base pore model */
		if ( sixMer.find('N') != std::string::npos ){
			generatedSignal.push_back( 68.0 );
		}
		else{
			generatedSignal.push_back( ( baseModel.at( sixMer ) ).first );
		}
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

	int numOfEvents = events.size();
	int numOfBases = generatedSignal.size();

	/*allocate the dynamic time warping lattice */
	std::vector< std::vector< double > > dtw( numOfEvents, std::vector< double >( numOfBases, std::numeric_limits< double >::max() ) );

	/*set the width of the Sakoe-Chiba band */
	int T = 25;

	/*initialisation */
	double runningSum = 0.0;
	for ( unsigned int i = 0; i < T; i++ ){

		runningSum += std::abs( events[ i ] - generatedSignal[ 0 ] );
		dtw[ i ][ 0 ] = runningSum;

	}

	runningSum = 0.0;
	for ( unsigned int i = 0; i < T; i++ ){

		runningSum += std::abs( events[ 0 ] - generatedSignal[ i ] );
		dtw[ 0 ][ i ] = runningSum; 

	}

	/*recursion: fill in the dynamic time warping lattice */
	for ( int j = 1; j < numOfBases; j++ ){
		
		int SCbandLower = std::max( 1, j - T );
		int SCbandUpper = std::min( numOfEvents, j + T );

		for ( int i = SCbandLower; i < SCbandUpper; i++ ){

			if ( i == 1 or j == 1 ){
				dtw[ i ][ j ] = std::abs( events[ i ] - generatedSignal[ j ] ) + std::min( dtw[ i - 1 ][ j ], std::min( dtw[ i ][ j - 1 ], dtw[ i - 1 ][ j - 1 ] ) );			
			}
			else{
				dtw[ i ][ j ] = std::abs( events[ i ] - generatedSignal[ j ] ) + std::min( dtw[ i - 1 ][ j - 1 ], std::abs( events[ i - 1 ] - generatedSignal[ j ] ) + std::min( dtw[ i - 2 ][ j - 1 ], std::min( dtw[ i - 1 ][ j ], std::abs( events[ i - 1 ] - generatedSignal[ j ] ) + dtw[ i - 2 ][ j ] ) ) );

			}

		}
	
	}

	/*termination: calculate the optimal warping path */
	int j = numOfBases - 1;
	int i = numOfEvents - 1;
	std::map< int, std::vector< int > > readPosToEventPos;

	while ( i > 2 and j > 2 ){

		readPosToEventPos[ j ].push_back( i );
		//std::cout << j << " " << i << std::endl;

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
			i--;
			j--;
		}
		else if ( m == 1 ){
			j--;
		}
		else if ( m == 2 ){
			i--;
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
