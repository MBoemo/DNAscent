//----------------------------------------------------------
// Copyright 2017 University of Oxford
// Written by Michael A. Boemo (michael.boemo@path.ox.ac.uk)
//----------------------------------------------------------

#ifndef DATA_STRUCTURES_H
#define DATA_STRUCTURES_H

#include <stdint.h>
#include <string.h> /*memset */
#include "data_structures.h"


/*struct for a matrix-like object that uses a contiguous block of memory */
template< typename type >
struct Matrix {
	type *cells;
	uint32_t numRows;
	uint32_t numCols;
};


/*for now, restrict attention to matrices of doubles and ints */
typedef Matrix< double > DoubleMatrix;
typedef Matrix< int > IntMatrix;


/*allocate a contiguous block of memory for the matrix */
template< typename type >
void allocate( Matrix< type >& matrix, uint32_t rows, uint32_t cols )
{
	matrix.numRows = rows;
	matrix.numRows = cols;
    
	uint32_t N = matrix.numRows * matrix.numCols;
	matrix.cells = ( type* ) malloc( N * sizeof( type ) );
	memset( matrix.cells, 0, N * sizeof( type ) );
}


/*return the appropriate cell, given row and column indicies */
template< typename type > 
inline uint32_t cell( const Matrix< type >& matrix, uint32_t row, uint32_t col )
{
	return row * matrix.numCols + col;
}


/*set a cell in the matrix to a value, using its row and column indices */
template< typename type, typename U >
inline void set( Matrix< type >& matrix, uint32_t row, uint32_t col, U value )
{
	uint32_t c = cell( matrix, row, col ); /*grab the cell that we're interested in */
	matrix.cells[ c ] = value; /*set that cell to the value that we want */
}


/*retrive a value from the matrix at a certain cell position, given by row and column indicies */
template< typename type >
inline type get( const Matrix< type >& matrix, uint32_t row, uint32_t col )
{
	uint32_t c = cell( matrix, row, col );
	return matrix.cells[ c ];
}

#endif
