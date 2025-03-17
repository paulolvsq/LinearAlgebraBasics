#include "LinearAlgebraBasics.h"

double *parallel_vector_matrix_product(double *A, int A_rows, int A_columns, double *X, int dimension) {

    if (A_rows <= 0 || A_columns <= 0 || dimension <= 0) {
        fprintf(stderr, "Error: Invalid dimensions for parallel vector-matrix product (A_rows=%d, A_columns=%d, dimension=%d). All dimensions must be strictly positive.\n",
                A_rows, A_columns, dimension);
        return NULL;
    }

    if (A_columns != dimension) {
        fprintf(stderr, "Error: Dimension mismatch. Columns of matrix (%d) must equal size of vector (%d).\n",
                A_columns, dimension);
        return NULL;
    }

    if (!A || !X) {
        fprintf(stderr, "Error: Null pointer detected for input matrix or vector in parallel_vector_matrix_product.\n");
        return NULL;
    }
    
    double *vector = malloc(A_rows * sizeof(double));

    if (!vector) {
        fprintf(stderr, "Error: Memory allocation failed for result vector in parallel_vector_matrix_product.\n");
        return NULL;
    }

    double value = 0.0;
#pragma omp parallel for
    for (int i = 0; i < A_rows; i++) {
	value = 0.0;
#pragma omp parallel for reduction(+:value)
	for (int j = 0; j < A_columns; j++) {
	    value += A[i * A_columns + j] * X[j];
	}
	vector[i] = value;
    }

    return vector;

}
