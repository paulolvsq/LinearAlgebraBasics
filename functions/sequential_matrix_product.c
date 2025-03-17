#include "LinearAlgebraBasics.h"

double *sequential_matrix_product(double *P, int P_rows, int P_columns, double *Q, int Q_rows, int Q_columns) {
    
    if (P_rows <= 0 || P_columns <= 0 || Q_rows <= 0 || Q_columns <= 0) {
        fprintf(stderr, "Error: Invalid dimensions provided. All dimensions must be strictly positive. "
                        "(P_rows=%d, P_columns=%d, Q_rows=%d, Q_columns=%d)\n", 
                        P_rows, P_columns, Q_rows, Q_columns);
        return NULL;
    }

    if (P_columns != Q_rows) {
        fprintf(stderr, "Error: Dimension mismatch. Columns of first matrix (%d) must equal rows of second matrix (%d).\n", 
                        P_columns, Q_rows);
        return NULL;
    }

    if (!P || !Q) {
        fprintf(stderr, "Error: Null pointer detected in sequential_matrix_product.\n");
        return NULL;
    }

    double *matrix = malloc(sizeof(double) * (P_rows * Q_columns));

    if (!matrix) {
        fprintf(stderr, "Error: Memory allocation failed for result matrix of size %dx%d.\n", 
                        P_rows, Q_columns);
        return NULL;
    }

    double value;
    
    for (int i = 0; i < P_rows; i++) {
	for (int j = 0; j < Q_columns; j++) {
	    value = 0;
	    for (int k = 0; k < P_columns; k++) {
		value += P[i * P_columns + k] * Q[k * Q_columns + j];
	    }
	    matrix[i * Q_columns + j] = value;
	}
    }

    return matrix;

}
