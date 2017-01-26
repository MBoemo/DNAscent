#include <iostream>
#include <vector>
#include <string>
#include <map>
#include <algorithm> /*count */
#include <utility> /*pair */
#include <iterator> /*distance */

/*modules */
#include "probability.h" /*Distribution, NormalDistribution, UniformDistribution */
#include "data_structures.h" /*DoubleMatrix, IntMatrix */

/*may want to rethink this later... but for now just use the preprocessor */
#define STARTNAME "START"
#define ENDNAME "END"

/*
##########################################################################################################
CLASS: State
EXAMPLES:
- state with a normal distribution: State(NormalDistribution(1,0), "unique_state_name1");
- state with a uniform distribution: State(UniformDistribution(0.5), "unique_state_name2",5.0);
- silent state: State(NULL, "unique_state_name3");
*/
class State {
	public:
		std::string name;
		double weight;
		Distribution *dist;
		State( Distribution *, std::string, double );
};

State::State( Distribution *x, std::string y, double z = 1.0 ){
	dist = x;
	name = y;
	weight = z;
}


/*
##########################################################################################################
CLASS: HiddenMarkovModel
EXAMPLE: HiddenMarkovModel( "AGTCG" );
*/
class HiddenMarkovModel {
	private:
		std::vector< State > startendStates;
		std::map< std::pair< std::string, std::string >, double > transitions;
		std::string reference;
		int silentStart;
		bool finalised;

	public:
		std::vector< State > states; /*TOCHANGE: For testing and development.  This should really be private. */
		std::vector< State > silentStates; /*TOCHANGE: For testing and development.  This should really be private. */
		HiddenMarkovModel( std::string );
		~HiddenMarkovModel( void );

		/*aux functions */
		void add_state( State );
		void add_transition( State, State, double );
		void finalise( );
		void sort_silentStates ( const int, const int );

		/*machine learning functions */
		std::pair< double, DoubleMatrix > forward( std::vector< double > );
		std::pair< double, DoubleMatrix > backward( std::vector< double > );
		std::vector< std::string > viterbi( std::vector< double > );
		void BaumWelch( std::vector< double >, int ); 
};

HiddenMarkovModel::HiddenMarkovModel( std::string x ) {
/*on constructor call, assign the reference, preallocate state vector memory, and set finalised boolean */

	reference = x;
    	/*pre-allocate enough memory for all states:
	  - each module contains six states: 6*x.length(),
	  - each thymidine requires two additional modules: 2*6*T_count
	  - start and end states: 2 */
    	int T_count = count( x.begin(), x.end(), 'T' );
	states.reserve( 6 * x.length() + 2 * 6 * T_count + 2 );

	/*pre-allocate memory for silent states - three (SS, D, SE) per module */
	silentStates.reserve( 3 * x.length() + 2 * 3 * T_count );
	
	/* add start and end states to the model */
	startendStates.reserve( 2 );
	State start( NULL, STARTNAME );
	State end( NULL, ENDNAME );
	startendStates.push_back( start );
	startendStates.push_back( end );

	/*on constructor call, initialise the finalised boolean to false to guard against calling algorithms before
	proper sorting has taken place */
	finalised = false;
}

HiddenMarkovModel::~HiddenMarkovModel( void ){
/*TOCHANGE: Destructor doesn't actually do anything yet. Fix later. */
}

void HiddenMarkovModel::add_state( State x ) {
/*add a state object to the states vector if it's emitting or silentStates if it's silent */
	if (x.dist){
	 	states.push_back( x );
	}
	else {
		silentStates.push_back( x );
	}
}

void HiddenMarkovModel::add_transition( State x, State y, double z ) {
/* keys a transition probability (double precision) of state x to state y with a pair of the two state names */
	std::pair< std::string, std::string > key = make_pair( x.name, y.name );
	transitions[key] = z;
}

void HiddenMarkovModel::finalise( ){
/* run after all states and transitions have been added, but before any HMM algorithms (forward, Viterbi, Baum-Welch) are executed */

	/*run a quicksort on the silent states so that all transitions between silent states go from low index to high index */
	HiddenMarkovModel::sort_silentStates(0, silentStates.size() - 1 );

	/*set the index in states where the silent states will begin */
	silentStart = states.size();

	/*append the (sorted) silent states vector to the end of the states vector */
	states.insert( states.end(), silentStates.begin(), silentStates.end() );

	/*flip the finalised boolean to true, where we can now run machine learning algorithms on the hmm object */
	finalised = true;
}

void HiddenMarkovModel::sort_silentStates( const int left, const int right ){
/*performs a quicksort on silent states so that transitions between silent states always to to a state with a higher index.  Simplifies forward, backward, and Viterbi algorithms. */

	std::pair< std::string, std::string > key;
	std::pair< std::string, std::string > rev_key;

	/*termination */
	if (left >= right) {
		return;
	}	

	/*partition */
	const int mid = left + (right - left) / 2;
   	const std::string pivot = silentStates[ mid ].name;
	
	/*swap the pivot to the left position */
 	std::swap( silentStates[ mid ], silentStates[ left ] );
  	int i = left + 1;
  	int j = right;
   	
	/*find two candidates for swapping, and swap them */
	while (i <= j) {
		
		rev_key = make_pair( pivot, silentStates[ i ].name );
        	while( i <= j && transitions.count( rev_key ) != 1 ){
        		i++;
			rev_key = make_pair( pivot, silentStates[ i ].name );
        	}

		key = make_pair( pivot, silentStates[ j ].name );
        	while(i <= j && transitions.count( key ) == 1 ) {
			j--;
			key = make_pair( pivot, silentStates[ j ].name );
        	}

        	if (i < j) {
            		std::swap( silentStates[ i ], silentStates[ j ] );
        	}
	}
	
	/*swap the pivot back to its original position */
	std::swap( silentStates[ i - 1 ], silentStates[ left ] );
	
	/*recursion */
	HiddenMarkovModel::sort_silentStates( left, i - 2 );
	HiddenMarkovModel::sort_silentStates( i, right );
}


std::pair< double, DoubleMatrix > HiddenMarkovModel::forward( std::vector <double> observations ) {

	/*guard - fall over if you call the forward algorithm before the hmm has been finalised */
	if (finalised == false){
		std::cout << "Exited with error. Finalise hidden Markov model object before running the forward algorithm." << std::endl;
		exit( -1 );
	}

	DoubleMatrix forward;
	allocate( forward, states.size(), observations.size() );
	std::pair< std::string, std::string > key;
	double logProbability, pdfEval;

	/*FUTURE: add something for moving between silent states before we see the first observation */

	/*INITIALISATION */
	/*emitting states */
	for ( int i = 0; i < silentStart; i++ ){

		key = make_pair( STARTNAME, states[ i ].name );

		if ( transitions.count( key ) > 0 ){
			logProbability = elnproduct( eln( transitions[ key ] ), eln( states[ i ].dist -> pdf( observations[ 0 ] ) ) );
		}
		else {
			logProbability = eln( 0.0 );
		}
		set( forward, i, 0, logProbability );
	}
	/*silent states */
	for ( int i = silentStart; i < states.size(); i++ ){

		key = make_pair( STARTNAME, states[ i ].name );

		if ( transitions.count( key ) > 0 ){
			logProbability = eln( transitions[ key ] );
		}
		else {
			logProbability = eln( 0.0 );
		}
		set( forward, i, 0, logProbability );
	}

	/*RECURSION (complexity is O(T*N^2) where T is the number of observations and N is the number of states) */
	for ( int t = 1; t < observations.size(); t++ ){
		
		/*handle emitting states */
		for ( int i = 0; i < silentStart; i++ ){ 
			
			pdfEval = states[ i ].dist -> pdf( observations[ t ] );
			logProbability = eln( 0.0 );
		
			for ( int k = 0; k < states.size(); k++ ){

				key = make_pair( states[ k ].name, states[ i ].name );
				
				if ( transitions.count( key ) > 0 ){
					logProbability = elnsum( logProbability, elnproduct( elnproduct( get( forward, k, t - 1 ), eln( transitions[ key ] ) ), eln( pdfEval ) ) );
				}
			}
			set( forward, i, t, logProbability );
		}
		/*handle silent states */
		for ( int i = silentStart; i < states.size(); i++ ){
			logProbability = eln( 0.0 );

			/*for all emitting states... */
			for ( int k = 0; k < silentStart; k++ ){

				key = make_pair( states[ k ].name, states[ i ].name );

				if ( transitions.count( key ) > 0 ){
					logProbability = elnsum( logProbability, elnproduct( get( forward, k, t ), eln( transitions[ key ] ) ) );
				}
			}
			/*for all silent states... */
			for ( int k = silentStart; k < i; k++ ){
			
				key = make_pair( states[ k ].name, states[ i ].name );

				if ( transitions.count( key ) > 0 ){
					logProbability = elnsum( logProbability, elnproduct( get( forward, k, t ), eln( transitions[ key ] ) ) );
				}
			}
			set( forward, i, t, logProbability );
		}
	}

	/*TERMINATION */
	double forwardProb = eln( 0.0 );
	for ( int i = 0; i < states.size(); i++ ){

		key = make_pair( states[ i ].name ,ENDNAME );	

		if ( transitions.count( key ) > 0) {
			forwardProb = elnsum( forwardProb, elnproduct( get( forward, i, observations.size() - 1), eln( transitions[ key ] ) ) );
		}
	}

	return std::make_pair( forwardProb, forward );
}


std::pair< double, DoubleMatrix > HiddenMarkovModel::backward( std::vector <double> observations ) {

	/*guard - fall over if you call the backward algorithm before the hmm has been finalised */
	if (finalised == false){
		std::cout << "Exited with error. Finalise hidden Markov model object before running the backward algorithm." << std::endl;
		exit( -1 );
	}

	DoubleMatrix backward;
	allocate( backward, states.size(), observations.size() );
	std::pair< std::string, std::string > key;
	double logProbability, pdfEval;

	/*FUTURE: add something for moving between silent states before we see the first observation */

	/*INITIALISATION */
	/*emitting and silent states (you can do these together since the end state doesn't have an emission) */
	for ( int i = 0; i < states.size(); i++ ){

		key = make_pair( states[ i ].name, ENDNAME );

		if ( transitions.count( key ) > 0 ){
			logProbability = eln( transitions[ key ] );
		}
		else {
			logProbability = eln( 0.0 );
		}
		set( backward, i, observations.size() - 1, logProbability);
	}

	/*RECURSION (complexity is O(T*N^2) where T is the number of observations and N is the number of states) */
	for (int t = observations.size() - 2; t >= 0; t--){
		
		/*handle emitting states */
		for ( int i = 0; i < silentStart; i++ ){ 
			
			pdfEval = states[ i ].dist -> pdf( observations[ t + 1 ] );
			logProbability = eln( 0.0 );
		
			for ( int k = 0; k < states.size(); k++ ){

				key = make_pair( states[ k ].name, states[ i ].name );
				
				if ( transitions.count( key ) > 0 ){
					logProbability = elnsum( logProbability, elnproduct( elnproduct( get( backward, k, t + 1 ), eln( transitions[ key ] ) ), eln( pdfEval ) ) );
				}
			}
			set( backward, i, t, logProbability );
		}
		/*handle silent states */
		for ( int i = silentStart; i < states.size(); i++ ){
			logProbability = eln( 0.0 );

			/*for all emitting states... */
			for ( int k = 0; k < silentStart; k++ ){

				key = make_pair( states[ k ].name, states[ i ].name );

				if ( transitions.count( key ) > 0 ){
					logProbability = elnsum( logProbability, elnproduct( get( backward, k, t ), eln( transitions[ key ] ) ) );
				}
			}
			/*for all silent states... */
			for ( int k = silentStart; k < i; k++ ){
			
				key = make_pair( states[ k ].name, states[ i ].name );

				if ( transitions.count( key ) > 0 ){
					logProbability = elnsum( logProbability, elnproduct( get( backward, k, t ), eln( transitions[ key ] ) ) );
				}
			}
			set( backward, i, t, logProbability );
		}
	}

	/*TERMINATION */
	double backwardProb = eln( 0.0 );
	for ( int i = 0; i < states.size(); i++ ){

		key = make_pair( STARTNAME, states[ i ].name );	

		if ( transitions.count( key ) > 0) {
			backwardProb = elnsum( backwardProb, elnproduct( get( backward, i, 0 ), eln( transitions[ key ] ) ) );
		}
	}

	return std::make_pair( backwardProb, backward );
}


std::vector< std::string > HiddenMarkovModel::viterbi( std::vector< double > observations ){

	DoubleMatrix viterbiLattice; /*main two dimensional lattice */
	allocate( viterbiLattice, states.size(), observations.size() );

	IntMatrix backtrace; /*stores state indices for the Viterbi backtrace */
	allocate( backtrace, states.size(), observations.size() );

	std::vector< double > maxViterbiLattice, maxBacktrace; /*helper vectors that we'll take the max element of */
	maxViterbiLattice.reserve( states.size() );
	maxBacktrace.reserve( states.size() );

	std::vector< std::string > viterbiFinal; /*final vector of the state names, in order, that most likely produced the observation */
	viterbiFinal.reserve( observations.size() );

	std::pair< std::string, std::string > key;
	double logProbability, pdfEval;
	int argMax;
	
	/*add something for moving between silent states before we see the first observation */

	/*INITIALISATION */
	/*emitting states */
	for ( int i = 0; i < silentStart; i++ ){

		key = make_pair( STARTNAME, states[i].name );

		if ( transitions.count( key ) > 0 ){
			logProbability = transitions[ key ] * states[i].dist -> pdf( observations[ 0 ] );
		}
		else {
			logProbability = 0;
		}

		set( viterbiLattice, i, 0, logProbability );
		set( backtrace, i, 0, 0);
	}
	/*silent states */
	for ( int i = silentStart; i < states.size(); i++ ){

		key = make_pair( STARTNAME, states[i].name );

		if ( transitions.count( key ) > 0 ){
			logProbability = transitions[ key ];
		}
		else {
			logProbability = 0;
		}

		set( viterbiLattice, i, 0, logProbability);
		set( backtrace, i, 0, 0);
	}

	/*RECURSION */
	for ( int t = 1; t < observations.size(); t++ ){

		/*emitting states */		
		for ( int i = 0; i < silentStart; i++ ){
		
			for ( int k = 0; k < states.size(); k++ ){
				
				key = make_pair( states[ k ].name, states[ i ].name );			
				maxViterbiLattice[ k ] = get( viterbiLattice, k, t - 1 ) * transitions[ key ] * states[ i ].dist -> pdf( observations[ t ] );
				maxBacktrace[ k ] = get( viterbiLattice, k, t - 1 ) * transitions[ key ];
			}

			/*the Viterbi lattice at (i,t) is the max element of the maxViterbiLattice vector */
			logProbability = *max_element( maxViterbiLattice.begin(), maxViterbiLattice.end() ); 
			set( viterbiLattice, i, t, logProbability);

			/*number of steps to get from begin to the max element, which is the state index (argmax) */
			argMax = distance( maxBacktrace.begin(), max_element( maxBacktrace.begin(), maxBacktrace.end() ) ); 
			set( backtrace, i, t, argMax );
		}
		/*silent states */
		for ( int i = silentStart; i < states.size(); i++ ){
			
			/*for all emitting states... */
			for ( int k = 0; k < silentStart; k++ ){
				
				key = make_pair( states[ k ].name, states[ i ].name );			
				maxViterbiLattice[ k ] = get( viterbiLattice, k, t ) * transitions[ key ];
				maxBacktrace[ k ] = get( viterbiLattice, k, t ) * transitions[ key ];
			}

			/*for all silent states... */
			for ( int k = silentStart; k < i; k++ ){
				
				key = make_pair( states[ k ].name, states[ i ].name );			
				maxViterbiLattice[ k ] = get( viterbiLattice, k, t ) * transitions[ key ];
				maxBacktrace[ k ] = get( viterbiLattice, k, t ) * transitions[ key ];
			}

			/*the Viterbi lattice at (i,t) is the max element of the maxViterbiLattice vector */
			logProbability = *max_element( maxViterbiLattice.begin(), maxViterbiLattice.end() ); 
			set( viterbiLattice, i, t, logProbability);

			/*number of steps to get from begin to the max element, which is the state index (argmax) */
			argMax = distance( maxBacktrace.begin(), max_element( maxBacktrace.begin(), maxBacktrace.end() ) ); 
			set( backtrace, i, t, argMax ); 
		}
	}

	/*TERMINATION */
	int tracebackIndex;
	for ( int i = 0; i < states.size(); i++ ){
			
		key = make_pair( states[ i ].name, ENDNAME );			
		maxViterbiLattice[ i ] = get( viterbiLattice, i, observations.size() ) * transitions[ key ]; /*maxBacktrace is equal to maxViterbiLattice in the termination step - just use this one for both */
	}

	tracebackIndex = distance( maxViterbiLattice.begin(), max_element( maxViterbiLattice.begin(), maxViterbiLattice.end() ) );

	/*TRACEBACK */
	/*traceback - initialise */
	viterbiFinal[ observations.size() - 1 ] = states[ tracebackIndex ].name;

	/*traceback - recursion */
	for (int t = observations.size() - 2; t > -1; t--){

		tracebackIndex = get( backtrace, tracebackIndex, t + 1 );
		viterbiFinal[ t ] = states[ tracebackIndex ].name;
	}

	return viterbiFinal;
}


void HiddenMarkovModel::BaumWelch( std::vector< double > observation, int maxIter ){
	/*gamma: probability of being in state j at time t */


	/*guard - fall over if you call the Baum-Welch training algorithm before the HMM has been finalised */
	if (finalised == false){
		std::cout << "Exited with error. Finalise hidden Markov model object before running the Baum-Welch algorithm." << std::endl;
		exit( -1 );
	}

	/*matrices for forward and backward algorithms */
	std::pair< double, DoubleMatrix > forwardPair = forward( observation );
	double forwardProb = forwardPair.first;
	DoubleMatrix alpha = forwardPair.second;
	std::pair< double, DoubleMatrix > backwardPair = backward( observation );
	DoubleMatrix beta = backwardPair.second;

	std::pair< std::string, std::string > key;
	double aNumerator, aDenominator, muNumerator, sigmaNumerator, bDenominator;

	/*main iteration loop */
	for ( int iter = 0; iter < maxIter; iter++ ){
	
		/*loop over states */
		for ( int j = 0; j < states.size(); j++ ){

			/*write something that only updates probability for normal distribution emitting states - not silent or uniform */

			/*EMISSIONS */
			/*for each state, update the mu, sigma values for each normal emitting state */
			muNumerator = 0.0;
			sigmaNumerator = 0.0;
			bDenominator = 0.0;
	
			for (int t = 0; t < observation.size(); t++){
				muNumerator = muNumerator + observation[ t ] * ( get( alpha, t, j) * get( beta, t, j) ) / forwardProb;
				sigmaNumerator = sigmaNumerator + pow( observation[ t ] - states[j].dist -> param1, 2 ) * ( get( alpha, t, j) * get( beta, t, j) ) / forwardProb;
				bDenominator = bDenominator + ( get( alpha, t, j) * get( beta, t, j) ) / forwardProb;
			}

			/*update the distribution parameters */
			states[j].dist -> param1 = muNumerator / bDenominator;
			states[j].dist -> param2 = sigmaNumerator / bDenominator;
			
			/*TRANSITIONS */
			/*update the transitions by looping over all states again, so we get all states x states */
			for ( int i = 0; i < states.size(); i++ ){

				aNumerator = 0.0;
				aDenominator = 0.0;

				key = make_pair( states[ i ].name, states[ j ].name );

				for ( int t = 0; t < observation.size() - 1; t++ ){


					aNumerator = aNumerator + ( get( alpha, t, i ) * transitions[ key ] * get( beta, t+1, j ) * states[j].dist -> pdf( observation[ t+1 ] ) ) / forwardProb;

					for ( int k = 0; k < states.size(); k++ ){
						aDenominator = aDenominator + ( get( alpha, t, i ) * transitions[key] * get( beta, t+1, k ) * states[k].dist -> pdf( observation[ t+1 ] ) ) / forwardProb;				}
				}
				
				transitions[ key ] = aNumerator / aDenominator; /*update transition sparse matrix */

			}

		}
	}
}
