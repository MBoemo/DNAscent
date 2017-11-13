//----------------------------------------------------------
// Copyright 2017 University of Oxford
// Written by Michael A. Boemo (michael.boemo@path.ox.ac.uk)
// This software is licensed under GPL-3.0.  You should have
// received a copy of the license with this software.  If
// not, please Email the author.
//----------------------------------------------------------


#include <algorithm>
#include <iostream>
#include <fstream>
#include "build_model.h"
#include "poreSpecificParameters.h"
#include "../Penthus/src/hmm.h"
#include "../Penthus/src/states.h"
#include "data_IO.h"


std::stringstream buildAndTrainHMM( std::string &reference, std::map< std::string, std::pair< double, double > > &basePoreModel, std::vector< std::vector< double > > &events, int &threads, bool verbose = false ){

	HiddenMarkovModel hmm = HiddenMarkovModel( 3*reference.length(), 3*reference.length() + 2 );

	std::pair< double, double > emissionMeanAndStd;

	/*STATES - vector (of vectors) to hold the states at each position on the reference - fill with dummy values */
	std::vector< std::vector< State > > states( 6, std::vector< State >( reference.length() - 5, State( NULL, "", "", "", 1.0 ) ) );

	/*DISTRIBUTIONS - vector to hold normal distributions, a single uniform and silent distribution to use for everything else */
	std::vector< NormalDistribution > nd;
	nd.reserve( reference.length() - 6 );		
	SilentDistribution sd( 0.0, 0.0 );
	UniformDistribution ud( 50.0, 130.0 );

	std::string loc;

	/*create the distributions that we need */			
	for ( unsigned int i = 0; i < reference.length() - 6; i++ ){

		emissionMeanAndStd = basePoreModel.at( reference.substr( i, 6 ) );		
		nd.push_back( NormalDistribution( emissionMeanAndStd.first, emissionMeanAndStd.second ) );		
	}

	/*add states to model, handle internal module transitions */
	for ( unsigned int i = 0; i < reference.length() - 6; i++ ){

		loc = std::to_string( i );
		std::string sixMer = reference.substr( i, 6 );

		states[ 0 ][ i ] = State( &sd, 		loc + "_SS",	sixMer,	"", 		1.0 );
		states[ 1 ][ i ] = State( &sd,		loc + "_D", 	sixMer,	"", 		1.0 );		
		states[ 2 ][ i ] = State( &ud,		loc + "_I", 	sixMer,	"", 		1.0 );
		states[ 3 ][ i ] = State( &nd[i], 	loc + "_M1", 	sixMer,	loc + "_match", 1.0 );
		states[ 4 ][ i ] = State( &nd[i], 	loc + "_M2", 	sixMer,	loc + "_match", 1.0 );
		states[ 5 ][ i ] = State( &sd, 		loc + "_SE", 	sixMer,	"", 		1.0 );

		/*add state to the model */
		for ( unsigned int j = 0; j < 6; j++ ){

			states[ j ][ i ].meta = sixMer;
			hmm.add_state( states[ j ][ i ] );
		}

		/*transitions between states, internal to a single base */
		/*from SS */
		hmm.add_transition( states[0][i], states[3][i], internalSS2M1 );
		hmm.add_transition( states[0][i], states[4][i], internalSS2M2 );

		/*from D */
		hmm.add_transition( states[1][i], states[2][i], internalD2I );

		/*from I */
		hmm.add_transition( states[2][i], states[2][i], internalI2I );
		hmm.add_transition( states[2][i], states[0][i], internalI2SS );

		/*from M1 */
		hmm.add_transition( states[3][i], states[3][i], internalM12M1 );
		hmm.add_transition( states[3][i], states[5][i], internalM12SE );

		/*from M2 */
		hmm.add_transition( states[4][i], states[4][i], internalM22M2 );
		hmm.add_transition( states[4][i], states[5][i], internalM22SE );

		/*from SE */
		hmm.add_transition( states[5][i], states[2][i], internalSE2I );		

	}

	/*add transitions between modules (external transitions) */
	for ( unsigned int i = 0; i < reference.length() - 7; i++ ){

		hmm.add_transition( states[1][i], states[1][i + 1], externalD2D );
		hmm.add_transition( states[1][i], states[0][i + 1], externalD2SS );
		hmm.add_transition( states[2][i], states[0][i + 1], externalI2SS );
		hmm.add_transition( states[5][i], states[0][i + 1], externalSE2SS );
		hmm.add_transition( states[5][i], states[1][i + 1], externalSE2D );
	}

	/*handle start states */
	hmm.add_transition( hmm.start, states[0][0], 0.5 );
	hmm.add_transition( hmm.start, states[1][0], 0.5 );

	/*handle end states */
	hmm.add_transition( states[1][reference.length() - 7], hmm.end, externalD2D + externalD2SS );
	hmm.add_transition( states[2][reference.length() - 7], hmm.end, externalI2SS );
	hmm.add_transition( states[5][reference.length() - 7], hmm.end, externalSE2SS + externalSE2D );

	hmm.finalise();

	hmm.BaumWelch( events, 1, 250, 0.0, false, threads, verbose );

	std::stringstream ss = hmm.summarise();

	return ss;
}


double buildAndDetectHMM( std::string &reference, std::map< std::string, std::pair< double, double > > &basePoreModel, std::map< std::string, std::pair< double, double > > &analogueModel, std::vector< double > &events, bool analogue ){

	HiddenMarkovModel hmm = HiddenMarkovModel( 3*reference.length() + 2, 3*reference.length() + 2 );

	std::pair< double, double > emissionMeanAndStd;

	/*STATES - vector (of vectors) to hold the states at each position on the reference - fill with dummy values */
	std::vector< std::vector< State > > states( 6, std::vector< State >( reference.length() - 5, State( NULL, "", "", "", 1.0 ) ) );

	/*DISTRIBUTIONS - vector to hold normal distributions, a single uniform and silent distribution to use for everything else */
	std::vector< NormalDistribution > nd;
	nd.reserve( reference.length() - 6 );
	SilentDistribution sd( 0.0, 0.0 );
	UniformDistribution ud( 50.0, 130.0 );

	/*garbage collection states to make up for inaccuracies in dynamic time warping */
	State gcStart( &ud, "gcStart", "", "", 1.0 );
	State gcEnd( &ud, "gcEnd", "", "", 1.0 );
	hmm.add_state( gcStart );
	hmm.add_state( gcEnd );

	std::string loc;

	/*create the distributions that we need */	
	for ( unsigned int i = 0; i < reference.length() - 6; i++ ){
		
		/*fill in analogue data in the appropriate position, if we have it.  otherwise, use the ONT pore model */
		if ( analogue == true and i == 7 and analogueModel.count( ( reference.substr( i, 6 ) ).replace( 3, 1, "B" ) ) > 0 ) {

			emissionMeanAndStd = analogueModel.at( ( reference.substr( i, 6 ) ).replace( 3, 1, "B" ) );

		}
		else if ( analogue == true and i == 8 and analogueModel.count( ( reference.substr( i, 6 ) ).replace( 2, 1, "B" ) ) > 0 ){

			emissionMeanAndStd = analogueModel.at( ( reference.substr( i, 6 ) ).replace( 2, 1, "B" ) );			

		}
		else{

			emissionMeanAndStd = basePoreModel.at( reference.substr( i, 6 ) );
		}

		nd.push_back( NormalDistribution( emissionMeanAndStd.first, emissionMeanAndStd.second ) );

	}

	/*add states to model, handle internal module transitions */
	for ( unsigned int i = 0; i < reference.length() - 6; i++ ){

		loc = std::to_string( i );

		states[ 0 ][ i ] = State( &sd, 		loc + "_SS", 	reference.substr( i, 6 ),	"", 		1.0 );
		states[ 1 ][ i ] = State( &sd,		loc + "_D", 	reference.substr( i, 6 ),	"", 		1.0 );		
		states[ 2 ][ i ] = State( &ud,	 	loc + "_I", 	reference.substr( i, 6 ),	"", 		1.0 );
		states[ 3 ][ i ] = State( &nd[ i ], 	loc + "_M1", 	reference.substr( i, 6 ),	loc + "_match", 1.0 );
		states[ 4 ][ i ] = State( &nd[ i ], 	loc + "_M2", 	reference.substr( i, 6 ),	loc + "_match", 1.0 );
		states[ 5 ][ i ] = State( &sd, 		loc + "_SE", 	reference.substr( i, 6 ),	"", 		1.0 );

		/*add state to the model */
		for ( unsigned int j = 0; j < 6; j++ ){

			states[ j ][ i ].meta = reference.substr( i, 6 );
			hmm.add_state( states[ j ][ i ] );
		}

		/*transitions between states, internal to a single base */
		/*from SS */
		hmm.add_transition( states[0][i], states[3][i], internalSS2M1 );
		hmm.add_transition( states[0][i], states[4][i], internalSS2M2 );

		/*from D */
		hmm.add_transition( states[1][i], states[2][i], internalD2I );

		/*from I */
		hmm.add_transition( states[2][i], states[2][i], internalI2I );
		hmm.add_transition( states[2][i], states[0][i], internalI2SS );

		/*from M1 */
		hmm.add_transition( states[3][i], states[3][i], internalM12M1 );
		hmm.add_transition( states[3][i], states[5][i], internalM12SE );

		/*from M2 */
		hmm.add_transition( states[4][i], states[4][i], internalM22M2 );
		hmm.add_transition( states[4][i], states[5][i], internalM22SE );

		/*from SE */
		hmm.add_transition( states[5][i], states[2][i], internalSE2I );		
	}

	/*add transitions between modules (external transitions) */
	for ( unsigned int i = 0; i < reference.length() - 7; i++ ){

		hmm.add_transition( states[ 1 ][ i ], states[ 1 ][ i + 1 ], externalD2D );
		hmm.add_transition( states[ 1 ][ i ], states[ 0 ][ i + 1 ], externalD2SS );
		hmm.add_transition( states[ 2 ][ i ], states[ 0 ][ i + 1 ], externalI2SS );
		hmm.add_transition( states[ 5 ][ i ], states[ 0 ][ i + 1 ], externalSE2SS );
		hmm.add_transition( states[ 5 ][ i ], states[ 1 ][ i + 1 ], externalSE2D );
	}

	/*handle start states */
	hmm.add_transition( hmm.start, gcStart, 1.0 );
	hmm.add_transition( gcStart, gcStart, 0.8 );
	hmm.add_transition( gcStart, states[ 0 ][ 0 ], 0.1 );
	hmm.add_transition( gcStart, states[ 1 ][ 0 ], 0.1 );

	/*handle end states */
	hmm.add_transition( states[ 1 ][ reference.length() - 7 ], gcEnd, externalD2D + externalD2SS );
	hmm.add_transition( states[ 2 ][ reference.length() - 7 ], gcEnd, externalI2SS );
	hmm.add_transition( states[ 5 ][ reference.length() - 7 ], gcEnd, externalSE2SS + externalSE2D );
	hmm.add_transition( gcEnd, gcEnd, 0.8 );
	hmm.add_transition( gcEnd, hmm.end, 0.2 );

	hmm.finalise();

	double seqProb = hmm.sequenceProbability( events );

	return seqProb;

}
