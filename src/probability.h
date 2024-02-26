//----------------------------------------------------------
// Copyright 2019 University of Oxford
// Written by Michael A. Boemo (mb915@cam.ac.uk)
// This software is licensed under GPL-3.0.  You should have
// received a copy of the license with this software.  If
// not, please Email the author.
//----------------------------------------------------------

#include <exception>

class NegativeLog : public std::exception {
	public:
		virtual const char * what () const throw () {
			return "Negative value passed to natural log function.";
		}
};


class DivideByZero : public std::exception {
	public:
		virtual const char * what () const throw () {
			return "lnQuot: Cannot divide by zero.";
		}
};

double eexp( double );
double eln( double );
double lnSum( double,  double );
double lnProd( double, double );
//double lnQuot( double, double );
bool lnGreaterThan( double, double );
double uniformPDF( double, double, double );
double normalPDF( double, double, double );
double cauchyPDF( double, double, double );
