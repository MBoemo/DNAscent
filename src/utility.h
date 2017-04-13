//----------------------------------------------------------
// Copyright 2017 University of Oxford
// Written by Michael A. Boemo (michael.boemo@path.ox.ac.uk)
//----------------------------------------------------------


/*utility functions for C++ Osiris */

std::string reverseComplement( std::string DNAseq ){

	std::reverse( DNAseq.begin(), DNAseq.end() );
	std::string revComp;

	for ( unsigned int i = 0; i < DNAseq.length(); i++ ){
	
		switch( DNAseq[ i ] ){
			case 'A' :
				revComp += 'T';
				break;
			case 'T' :
				revComp += 'A';
				break;
			case 'G' :
				revComp += 'C';
				break;
			case 'C' :
				revComp += 'G';
				break;
			default:
				std::cout << "Exiting with error.  Invalid character passed to reverse complement function.  Must be A, T, G, or C." << std::endl;
				exit( EXIT_FAILURE );
		}

	}

	return revComp;

}
