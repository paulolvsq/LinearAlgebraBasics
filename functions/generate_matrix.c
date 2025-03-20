#include "LinearAlgebraBasics.h"
#define DOUBLE 100.0

/**
 * @brief Generates a matrix with random double values between 0 and DOUBLE.
 *
 * This function creates a matrix of size rows x columns, where each element
 * is a random double value in the range [0, DOUBLE].
 *
 * @param rows Number of rows in the matrix (must be positive).
 * @param columns Number of columns in the matrix (must be positive).
 *
 * @return Pointer to the generated matrix on success, or NULL on failure due to invalid dimensions
 *         or memory allocation errors.
 *
 */

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

/**
 * @brief Generates an identity matrix of given dimension.
 *
 * This function creates a square identity matrix of size dimension x dimension,
 * where all diagonal elements are 1.0 and all off-diagonal elements are 0.0.
 *
 * @param dimension Dimension of the identity matrix (must be positive).
 *
 * @return Pointer to the generated identity matrix on success, or NULL on failure due to invalid dimension
 *         or memory allocation errors.
 */

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
	    
