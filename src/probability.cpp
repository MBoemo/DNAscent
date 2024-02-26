//----------------------------------------------------------
// Copyright 2019 University of Oxford
// Written by Michael A. Boemo (mb915@cam.ac.uk)
// This software is licensed under GPL-3.0.  You should have
// received a copy of the license with this software.  If
// not, please Email the author.
//----------------------------------------------------------

#include <algorithm>
#include <exception>
#include <vector>
#include "probability.h"
#include "error_handling.h"
#define _USE_MATH_DEFINES


/*
Functions for computing probabilities in log space to improve numerical stability.  Most algorithms
inspired by Numerically Stable Hidden Markov Model Implementation by Tobias P. Mann
(see http://bozeman.genome.washington.edu/compbio/mbt599_2006/hmm_scaling_revised.pdf).
*/

double eexp( double x ){
/*Map from log space back to linear space. */

	if ( std::isnan( x ) ) {
		return 0.0;
	}
	else {
		return exp( x );
	}
}


double eln( double x ){
/*Map from linear space to log space. */

	if (x == 0.0){
		return NAN;
	}
	else if (x > 0.0){
		return log( x );
	}
	else{
		throw NegativeLog();
	}
}


double lnSum( double ln_x, double ln_y ){
/*evalutes the quotient ln_x + ln_y */

	/*if one of the arguments is NAN, go into special handling for that */
	if ( std::isnan( ln_x ) || std::isnan( ln_y ) ){

        	if ( std::isnan( ln_x ) && std::isnan( ln_y ) ){
			return NAN;
		}
		else if ( std::isnan( ln_x ) ){
			return ln_y;
		}
		else{
			return ln_x;
		}
	}
	/*Otherwise, compute the output.  For numerical stability, we always want to subtract the larger argument from the smaller argument in the exponential term. */
	else{

		if ( ln_x > ln_y ){
			return ln_x + eln( 1.0 + eexp( ln_y - ln_x ) );
		}
		else{
			return ln_y + eln( 1.0 + eexp( ln_x - ln_y ) );
		}
	}
}


double lnProd( double ln_x, double ln_y ){
/*evalutes the quotient ln_x*ln_y */

	if ( std::isnan( ln_x ) || std::isnan( ln_y ) ){
		return NAN;
	}
	else{
		return ln_x + ln_y;
	}
}

/*
double lnQuot( double ln_x, double ln_y ){
//evalutes the quotient ln_x/ln_y

	if ( std::isnan( ln_y ) ){
        	throw DivideByZero();
	}
	else if ( std::isnan( ln_x ) ){
		return NAN;
	}
	else{
		return ln_x - ln_y;
	}
}
*/


bool lnGreaterThan( double ln_x, double ln_y ){
/*evalutes whether ln_x is greater than ln_y, and returns a boolean */

	if ( std::isnan( ln_x ) || std::isnan( ln_y ) ){

		if ( std::isnan( ln_x ) || std::isnan( ln_y ) == false ){
			return false;
		}
        	else if ( std::isnan( ln_x ) == false || std::isnan( ln_y ) ){
			return true;
		}
		else{
			return false;
		}
	}
	else{
		if ( ln_x > ln_y ){
			return true;
		}
		else{
			return false;
		}
	}
}


double uniformPDF( double lb, double ub, double x ){

	if ( x >= lb && x <= ub ){

		return 1.0/( ub - lb );
	}
	else {
		return 0.0;
	}
};


double normalPDF( double mu, double sigma, double x ){

	return ( 1.0/sqrt( 2.0*pow( sigma, 2.0 )*M_PI ) )*exp( -pow( x - mu , 2.0 )/( 2.0*pow( sigma, 2.0 ) ) );
}


double cauchyPDF( double loc, double scale, double x ){

	return 1./( (scale*M_PI) * ( 1. + pow((x-loc)/scale, 2.) ) );
}
