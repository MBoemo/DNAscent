//----------------------------------------------------------
// Copyright 2017 University of Oxford
// Written by Michael A. Boemo (michael.boemo@path.ox.ac.uk)
//----------------------------------------------------------

#include <iostream>
#include <math.h> /*pow, sqrt, isnan, NAN */

#define _USE_MATH_DEFINES

/*
##########################################################################################################
CLASS: LogSpace
Functions for computing probabilities in log space to improve numerical stability.  Most algorithms 
inspired by Numerically Stable Hidden Markov Model Implementation by Tobias P. Mann
(see http://bozeman.genome.washington.edu/compbio/mbt599_2006/hmm_scaling_revised.pdf).
*/

double eexp( double x ){
/* if we take the log of probability 0, then we get NAN.  But we should return 0 when we take the exponential to map it back from log space */
	if ( isnan( x ) ) {
		return 0.0;
	}	
	else {
		return exp( x );
	}
}

double eln( double x ){
	if (x == 0.0){
		return NAN;
	}
	else if (x > 0.0){
		return log( x );
	}	
	else{
		std::cout << "Exited with error.  Negative value passed to natural log function." << std::endl;
		exit( -1 );
	}
}

double elnsum( double ln_x, double ln_y ){
	/*TOCHANGE: This could be cleaned up, and I'm not sure about the handling when both ln_x and ln_y are NAN */
	if ( isnan( ln_x ) ){
		if ( isnan( ln_x ) ){
			return ln_y;
		}
		else{
			return ln_x;
		}
	}
	else{
		if ( ln_x > ln_y ){
			return ln_x + eln( 1.0 + eexp( ln_y - ln_x ) );
		}
		else{
			return ln_y + eln( 1.0 + eexp( ln_x - ln_y ) );
		}
	}
}

double elnproduct( double ln_x, double ln_y ){

	if ( isnan( ln_x ) || isnan( ln_y ) ){
		return NAN;
	}
	else{
		return ln_x + ln_y;
	}
}


/*
##########################################################################################################
BASE CLASS: Distribution
*/
class Distribution{
	protected:
		Distribution( double, double);

	public:
		double param1;
		double param2;
		virtual double pdf( double ) = 0;

};

Distribution::Distribution( double x, double y ){
	param1 = x;
	param2 = y;
}

/*
##########################################################################################################
DERIVED CLASS: NormalDistribution (BASE CLASS: Distribution)
EXAMPLE: NormalDistribution( 0, 1 );
*/
class NormalDistribution: public Distribution {
	

	public:
		NormalDistribution( double, double );
		double mu, sigma;
		virtual double pdf( double x ){
			return ( 1.0/sqrt( 2.0*pow( sigma, 2.0 )*M_PI ) )*exp( -pow( x - mu , 2.0 )/( 2.0*pow( sigma, 2.0 ) ) );
		} ;
};

NormalDistribution::NormalDistribution( double x, double y ): Distribution( x, y ) {
	mu = x;
	sigma = y;
}

/*
##########################################################################################################
DERIVED CLASS: UniformDistribution (BASE CLASS: Distribution)
EXAMPLE: UniformDistribution( 0.0, 5.0 );
*/
class UniformDistribution: public Distribution{
	double a, b;

	public:
		UniformDistribution( double, double );
		virtual double pdf( double x){
			double probability;	
	
			if (x >= a && x <= b){
				probability = 1.0/(b - a);
			}
			else {
				probability = 0.0;
			}

			return probability;
		}
};

UniformDistribution::UniformDistribution( double x, double y ): Distribution(x,y) {
	a = x;
	b = y;	
}

