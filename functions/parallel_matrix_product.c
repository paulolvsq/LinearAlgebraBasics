#include "LinearAlgebraBasics.h"

/**
 * @brief Computes the product of two matrices P and Q in parallel using OpenMP.
 *
 * This function calculates the matrix product C = P * Q using parallelization
 * to improve performance. Each element of the resulting matrix is computed
 * by summing the products of corresponding elements from rows of P and columns of Q.
 *
 * @param P Pointer to the first input matrix (size: P_rows x P_columns).
 * @param P_rows Number of rows in matrix P (must be positive).
 * @param P_columns Number of columns in matrix P (must equal Q_rows).
 * @param Q Pointer to the second input matrix (size: Q_rows x Q_columns).
 * @param Q_rows Number of rows in matrix Q (must equal P_columns).
 * @param Q_columns Number of columns in matrix Q (must be positive).
 *
 * @return Pointer to the resulting matrix (size: P_rows x Q_columns) on success,
 *         or NULL on failure due to invalid dimensions, null pointers, or memory allocation errors.
 */

double *parallel_matrix_product(double *P, int P_rows, int P_columns, double *Q, int Q_rows, int Q_columns) {

    if (P_rows <= 0 || P_columns <= 0 || Q_rows <= 0 || Q_columns <= 0) {
        fprintf(stderr, "Error: Invalid dimensions for parallel matrix product (P_rows=%d, P_columns=%d, Q_rows=%d, Q_columns=%d). All dimensions must be strictly positive.\n",
                P_rows, P_columns, Q_rows, Q_columns);
        return NULL;
    }

    if (P_columns != Q_rows) {
        fprintf(stderr, "Error: Dimension mismatch. Columns of first matrix (%d) must equal rows of second matrix (%d).\n",
                P_columns, Q_rows);
        return NULL;
    }

    if (!P || !Q) {
        fprintf(stderr, "Error: Null pointer detected for input matrices in parallel_matrix_product.\n");
        return NULL;
    }

    double *matrix = malloc(sizeof(double) * (P_rows * Q_columns));

    if (!matrix) {
        fprintf(stderr, "Error: Memory allocation failed for result matrix in parallel_matrix_product.\n");
        return NULL;
    }

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
