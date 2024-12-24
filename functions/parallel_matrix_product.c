#include "LinearAlgebraBasics.h"

double *parallel_matrix_product(double *P, int P_rows, int P_columns, double *Q, int Q_rows, int Q_columns) {
    
    // ASSERTION : THE USER WILL ALWAYS GIVE MATRICES OF THE RIGHT SIZE AS INPUT

    double *matrix = malloc(sizeof(double) * (P_rows * Q_columns));

    double value;
#pragma omp parallel for
    for (int i = 0; i < P_rows; i++) {
	for (int j = 0; j < Q_columns; j++) {
	    value = 0;
#pragma omp parallel for reduction(+:value)
	    for (int k = 0; k < P_columns; k++) {
		value += P[i * P_columns + k] * Q[k * Q_columns + j];
	    }
	    matrix[i * Q_columns + j] = value;
	}
    }

    return matrix;

}
