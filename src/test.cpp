//----------------------------------------------------------
// Copyright 2017 University of Oxford
// Written by Michael A. Boemo (michael.boemo@path.ox.ac.uk)
//----------------------------------------------------------

#include <iostream>

/* local header files */
#include "hmm.h"


int main(){

	/* TEST FOR ADDING START/END STATES 
	HiddenMarkovModel hmm("AGTC");
	for (int i = 0; i < hmm.startEnd.size(); i++){
	std::cout << hmm.startEnd[i].name << std::endl;
	}
	*/
	
	/* TEST FOR SORTING SILENT STATES */
	HiddenMarkovModel hmm("AGTC");
	hmm.add_state(State(NULL,"first_added",1.0));
	hmm.add_state(State(NULL,"second_added",1.0));
	hmm.add_state(State(NULL,"third_added",1.0));
	hmm.add_state(State(NULL,"fourth_added",1.0));
	hmm.add_transition(State(NULL,"second_added",1.0), State(NULL,"first_added",1.0), 1.0);
	hmm.add_transition(State(NULL,"third_added",1.0), State(NULL,"second_added",1.0), 1.0);
	hmm.add_transition(State(NULL,"fourth_added",1.0), State(NULL,"third_added",1.0), 1.0);
	for (int i = 0; i < hmm.silentStates.size(); i++){
		std::cout << hmm.silentStates[i].name << std::endl;
	}
	hmm.finalise();
	for (int i = 0; i < hmm.states.size(); i++){
		std::cout << hmm.states[i].name << std::endl;
	}


	/* TEST FOR PDF POLYMORPHISM 
	NormalDistribution nd(0.0,1.0);
	NormalDistribution *np = &nd;
	State state(np,"name",1.0);
	std::cout << state.dist -> pdf(0.0) << std::endl; 
	*/

	
	return 0;
}

