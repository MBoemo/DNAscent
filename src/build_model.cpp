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
	nd.reserve( reference.length() - 5 );		
	SilentDistribution sd( 0.0, 0.0 );
	UniformDistribution ud( 50.0, 130.0 );

	std::string loc;

	/*create the distributions that we need */			
	for ( unsigned int i = 0; i < reference.length() - 5; i++ ){

		emissionMeanAndStd = basePoreModel.at( reference.substr( i, 5 ) );		
		nd.push_back( NormalDistribution( emissionMeanAndStd.first, emissionMeanAndStd.second ) );		
	}

	/*add states to model, handle internal module transitions */
	for ( unsigned int i = 0; i < reference.length() - 5; i++ ){

		loc = std::to_string( i );
		std::string sixMer = reference.substr( i, 5 );

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
	for ( unsigned int i = 0; i < reference.length() - 6; i++ ){

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
	hmm.add_transition( states[1][reference.length() - 6], hmm.end, externalD2D + externalD2SS );
	hmm.add_transition( states[2][reference.length() - 6], hmm.end, externalI2SS );
	hmm.add_transition( states[5][reference.length() - 6], hmm.end, externalSE2SS + externalSE2D );

	hmm.finalise();

	hmm.BaumWelch( events, 1.0, 250, 0.0, false, threads, verbose );

	std::stringstream ss = hmm.summarise();

	return ss;
}


double buildAndDetectHMM( std::string &reference, std::map< std::string, std::pair< double, double > > &basePoreModel, std::map< std::string, std::pair< double, double > > &analogueModel, std::vector< double > &events, unsigned int windowLength, bool analogue ){

	HiddenMarkovModel hmm = HiddenMarkovModel( 3*reference.length(), 3*reference.length() + 2 );

	std::pair< double, double > emissionMeanAndStd;

	/*STATES - vector (of vectors) to hold the states at each position on the reference - fill with dummy values */
	std::vector< std::vector< State > > states( 6, std::vector< State >( reference.length() - 5, State( NULL, "", "", "", 1.0 ) ) );

	/*DISTRIBUTIONS - vector to hold normal distributions, a single uniform and silent distribution to use for everything else */
	std::vector< NormalDistribution > nd;
	nd.reserve( reference.length() - 5 );		
	SilentDistribution sd( 0.0, 0.0 );
	UniformDistribution ud( 50.0, 130.0 );

	std::string loc;
	std::vector< unsigned int > analoguePositions = {windowLength, windowLength-1, windowLength-2, windowLength-3}; //for positions 1,2,3,4 of a 5mer

	/*create the distributions that we need */	
	for ( unsigned int i = 0; i < reference.length() - 5; i++ ){

		/*use the trained analogue model if we're at the analogue position */
		if ( analogue and std::find(analoguePositions.begin(), analoguePositions.end(), i) != analoguePositions.end() ){

			emissionMeanAndStd = analogueModel.at((reference.substr(i, 5)).replace(windowLength - i, 1, "B"));
		}
		/*otherwise use the ONT model */
		else{

			emissionMeanAndStd = basePoreModel.at( reference.substr( i, 5 ) );
		}
		nd.push_back( NormalDistribution( emissionMeanAndStd.first, emissionMeanAndStd.second ) );

	}

	/*add states to model, handle internal module transitions */
	for ( unsigned int i = 0; i < reference.length() - 5; i++ ){

		loc = std::to_string( i );
		std::string fiveMer = reference.substr( i, 5 );

		states[ 0 ][ i ] = State( &sd, 		loc + "_SS",	fiveMer,	"", 		1.0 );
		states[ 1 ][ i ] = State( &sd,		loc + "_D", 	fiveMer,	"", 		1.0 );		
		states[ 2 ][ i ] = State( &ud,		loc + "_I", 	fiveMer,	"", 		1.0 );
		states[ 3 ][ i ] = State( &nd[i], 	loc + "_M1", 	fiveMer,	loc + "_match", 1.0 );
		states[ 4 ][ i ] = State( &nd[i], 	loc + "_M2", 	fiveMer,	loc + "_match", 1.0 );
		states[ 5 ][ i ] = State( &sd, 		loc + "_SE", 	fiveMer,	"", 		1.0 );

		/*add state to the model */
		for ( unsigned int j = 0; j < 6; j++ ){

			states[ j ][ i ].meta = fiveMer;
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
	for ( unsigned int i = 0; i < reference.length() - 6; i++ ){

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
	hmm.add_transition( states[1][reference.length() - 6], hmm.end, externalD2D + externalD2SS );
	hmm.add_transition( states[2][reference.length() - 6], hmm.end, externalI2SS );
	hmm.add_transition( states[5][reference.length() - 6], hmm.end, externalSE2SS + externalSE2D );

	hmm.finalise();

	double seqProb = hmm.sequenceProbability( events );

	return seqProb;

}
