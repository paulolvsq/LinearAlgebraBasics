#include "LinearAlgebraBasics.h"

double *sequential_vector_matrix_product(double *A, int A_rows, int A_columns, double *X, int dimension) {

    if (A_rows <= 0 || A_columns <= 0 || dimension <= 0) {
        fprintf(stderr, "Error: Invalid dimensions provided. All dimensions must be strictly positive. (A_rows=%d, A_columns=%d, dimension=%d)\n", 
                A_rows, A_columns, dimension);
        return NULL;
    }

    if (A_columns != dimension) {
        fprintf(stderr, "Error: Dimension mismatch. Matrix columns (%d) must equal vector size (%d).\n", 
                A_columns, dimension);
        return NULL;
    }

    if (!A || !X) {
        fprintf(stderr, "Error: Null pointer detected in sequential_vector_matrix_product.\n");
        return NULL;
    }

    double *vector = malloc(A_rows * sizeof(double));

    if (!vector) {
        fprintf(stderr, "Error: Memory allocation failed for result vector of size %d.\n", A_rows);
        return NULL;
    }

    double value = 0.0;
    
    for (int i = 0; i < A_rows; i++) {
	value = 0.0;
	for (int j = 0; j < A_columns; j++) {
	    value += A[i * A_columns + j] * X[j];
	}
	vector[i] = value;
    }

    return vector;

}
