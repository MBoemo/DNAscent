//----------------------------------------------------------
// Copyright 2019-2020 University of Oxford
// Written by Michael A. Boemo (mb915@cam.ac.uk)
// This software is licensed under GPL-3.0.  You should have
// received a copy of the license with this software.  If
// not, please Email the author.
//----------------------------------------------------------

#ifndef REGIONS_H
#define REGIONS_H

#include <cassert>

/*function prototypes */
int regions_main( int argc, char** argv );


class AnalogueScore{

	private:
		double _score = 0.0;
		bool _isSet = false;
	public:
		void set(double s){
			_score = s;
			_isSet = true;
		}
		double get(void){
			assert(_isSet);
			return _score;
		}
};

class Line{

	private:
		AnalogueScore BrdU;
		AnalogueScore Methyl;
		AnalogueScore BrdUvMethyl;
		unsigned int position;

	public:
		Line(unsigned int position){
			this -> position = position;
		}
		void setBrdUScore(double s){
			BrdU.set(s);
		}
		void setBrdUvMethylScore(double s){
			BrdUvMethyl.set(s);
		}
		void setMethylScore(double s){
			Methyl.set(s);
		}
		double getBrdUScore(void){
			return BrdU.get();
		}
		double getBrdUvMethylScore(void){
			return BrdUvMethyl.get();
		}
		double getMethylScore(void){
			return Methyl.get();
		}
};

#endif
