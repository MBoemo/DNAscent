//----------------------------------------------------------
// Copyright 2017 University of Oxford
// Written by Michael A. Boemo (michael.boemo@path.ox.ac.uk)
//----------------------------------------------------------


#include "build_model.h"

HiddenMarkovModel build_trainingHMM( std::string reference, std::map< std::string, std::pair< double, double > > basePoreModel ){

	/*preallocate states - we need len(ref)*3 emitting states and len(ref)*3 silent state s */	
	HiddenMarkovModel hmm( 3*reference.length(), 3*reference.length() + 2 );


	return hmm;


}
