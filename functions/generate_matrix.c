#include "LinearAlgebraBasics.h"
#define DOUBLE 100.0

double *generate_matrix_double (int rows, int columns) {

    if (rows <= 0 || columns <= 0) {
        fprintf(stderr, "Error: invalid dimensions (rows=%d, columns=%d). Both must be positive..\n", rows, columns);
        return NULL;
    }

    double *matrix = malloc((rows * columns) * sizeof(double));

    if (!matrix) {
        fprintf(stderr, "Error: memory allocation failed for a matrix of size %dx%d.\n", rows, columns);
        return NULL;
    }
    
    srand(time(NULL));

    for (int i = 0; i < rows; i++) 
	for (int j = 0; j < columns; j++) 
	    matrix[i * columns + j] = (double) (rand()) / RAND_MAX * DOUBLE;    

     return matrix;

}

double *generate_identity_matrix (int dimension) {

    if (dimension <= 0) {
        fprintf(stderr, "Error: Invalid dimension (%d). Must be strictly positive.\n", dimension);
        return NULL;
    }

    double *matrix = malloc((dimension * dimension) * sizeof(double));

    if (!matrix) {
        fprintf(stderr, "Error: Memory allocation failed for identity matrix of size %dx%d.\n", dimension, dimension);
        return NULL;
    }

    for (int i = 0; i < dimension; i++) {
	for (int j = 0; j < dimension; j++) {
	    if (i == j) matrix[i * dimension + j] = 1.0;
	    else matrix[i * dimension + j] = 0.0;
	}
    }
    
    return matrix;
    
}
	    
