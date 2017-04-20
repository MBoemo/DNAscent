//----------------------------------------------------------
// Copyright 2017 University of Oxford
// Written by Michael A. Boemo (michael.boemo@path.ox.ac.uk)
//----------------------------------------------------------


#include "build_model.h"
#include "poreSpecificParameters.h"


HiddenMarkovModel *build_trainingHMM( std::string &reference, std::map< std::string, std::pair< double, double > > &basePoreModel ){

	std::cout << "Building HMM... " << std::endl;

	HiddenMarkovModel *hmm = new HiddenMarkovModel( 3*reference.length(), 3*reference.length() + 2 );

	std::pair< double, double > emissionMeanAndStd;

	/*STATES - vector (of vectors) to hold the states at each position on the reference - fill with dummy values */
	static std::vector< std::vector< State > > states( 6, std::vector< State >( reference.length() - 5, State( NULL,"","", 1.0 ) ) );

	/*DISTRIBUTIONS - vector to hold normal distributions, a single uniform and silent distribution to use for everything else */
	static std::vector< NormalDistribution > nd;
	static SilentDistribution sd( 0.0, 0.0 );
	static UniformDistribution ud( 30.0, 130.0 );

	std::string loc;

	/*create the distributions that we need */	
	for ( unsigned int i = 0; i < reference.length() - 6; i++ ){
		
		emissionMeanAndStd = basePoreModel.at( reference.substr( i, 6 ) );
		nd.push_back( NormalDistribution( emissionMeanAndStd.first, emissionMeanAndStd.second ) );

	}

	/*add states to model, handle internal module transitions */
	for ( unsigned int i = 0; i < reference.length() - 6; i++ ){

		loc = std::to_string( i );

		states[ 0 ][ i ] = State( &sd, 		loc + "_SS", 	"", 		1.0 );
		states[ 1 ][ i ] = State( &sd,		loc + "_D", 	"", 		1.0 );		
		states[ 2 ][ i ] = State( &ud,	 	loc + "_I", 	"", 		1.0 );
		states[ 3 ][ i ] = State( &nd[ i ], 	loc + "_M1", 	loc + "_match", 1.0 );
		states[ 4 ][ i ] = State( &nd[ i ], 	loc + "_M2", 	loc + "_match", 1.0 );
		states[ 5 ][ i ] = State( &sd, 		loc + "_SE", 	"", 		1.0 );

		/*add state to the model */
		for ( unsigned int j = 0; j < 6; j++ ){
			hmm -> add_state( states[ j ][ i ] );
		}

		/*transitions between states, internal to a single base */
		/*from SS */
		hmm -> add_transition( states[0][i], states[3][i], internalSS2M1 );
		hmm -> add_transition( states[0][i], states[4][i], internalSS2M2 );

		/*from D */
		hmm -> add_transition( states[1][i], states[2][i], internalD2I );

		/*from I */
		hmm -> add_transition( states[2][i], states[2][i], internalI2I );
		hmm -> add_transition( states[2][i], states[0][i], internalI2SS );

		/*from M1 */
		hmm -> add_transition( states[3][i], states[3][i], internalM12M1 );
		hmm -> add_transition( states[3][i], states[5][i], internalM12SE );

		/*from M2 */
		hmm -> add_transition( states[4][i], states[4][i], internalM22M2 );
		hmm -> add_transition( states[4][i], states[5][i], internalM22SE );

		/*from SE */
		hmm -> add_transition( states[5][i], states[2][i], internalSE2I );		

	}

	/*add transitions between modules (external transitions) */
	for ( unsigned int i = 0; i < reference.length() - 7; i++ ){

		hmm -> add_transition( states[ 1 ][ i ], states[ 1 ][ i + 1 ], externalD2D );
		hmm -> add_transition( states[ 1 ][ i ], states[ 0 ][ i + 1 ], externalD2SS );
		hmm -> add_transition( states[ 2 ][ i ], states[ 0 ][ i + 1 ], externalI2SS );
		hmm -> add_transition( states[ 5 ][ i ], states[ 0 ][ i + 1 ], externalSE2SS );
		hmm -> add_transition( states[ 5 ][ i ], states[ 1 ][ i + 1 ], externalSE2D );

	}

	/*handle start states */
	hmm -> add_transition( hmm -> start, states[ 0 ][ 0 ], 0.5 );
	hmm -> add_transition( hmm -> start, states[ 1 ][ 0 ], 0.5 );

	/*handle end states */
	hmm -> add_transition( states[ 1 ][ reference.length() - 7 ], hmm -> end, externalD2D + externalD2SS );
	hmm -> add_transition( states[ 2 ][ reference.length() - 7 ], hmm -> end, externalI2SS );
	hmm -> add_transition( states[ 5 ][ reference.length() - 7 ], hmm -> end, externalSE2SS + externalSE2D );

	std::cout << "\tFinalising and sorting states..." << std::endl;
	hmm -> finalise();
	std::cout << "\tDone." << std::endl;

	std::cout << "Done." << std::endl;

	return hmm;

}
