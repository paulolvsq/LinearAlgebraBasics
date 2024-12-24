#include "LinearAlgebraBasics.h"

double *parallel_vector_matrix_product(double *A, int A_rows, int A_columns, double *X, double dimension) {

    double *vector = malloc(A_rows * sizeof(double));

    // ASSERTION : THE USER WILL ALWAYS GIVE MATRICES AND VECOTR OF THE RIGHT SIZE AS INPUT

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
