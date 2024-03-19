#include "LinearAlgebraBasics.h"

// ASSERTION : THE USER WILL ALWAYS GIVE MATRICES OF THE RIGHT SIZE AS INPUT 

double *matrices_addition(double *A, double *B, int rows, int columns) {

    double *matrix = malloc(rows * columns * sizeof(double));

    for (int i = 0; i < rows; i++)
	for (int j = 0; j < columns; j++)
	    matrix[i * columns + j] = A[i * columns + j] + B[i * columns + j];

    return matrix;

}

double *matrix_scalar_multiplication(double *A, int rows, int columns, int scalar) {

    double *matrix = malloc(rows * columns * sizeof(double));

    for (int i = 0; i < rows; i++)
	for (int j = 0; j < columns; j++)
	    matrix[i * columns + j] = scalar * A[i * columns + j];

    return matrix;
	       
}

double *matrix_transpose(double *A, int rows, int columns) {

    double *matrix = malloc(rows * columns * sizeof(double));

    for (int i = 0; i < columns; i++)
	for (int j = 0; j < rows; j++)
	    matrix[i * rows + j] = A[i * rows + j];

    return matrix;

}

double matrix_trace(double *A, int dimension) {

    double trace = 0.0;

    for (int i = 0; i < dimension; i++) 
	for (int j = 0; j < dimension; j++) 
	    if (i == j) trace += A[i * dimension + j];
	

    return trace;
    
}

double matrix_norm(double *A, int rows, int columns) {

    // 1-NORM

    double max = -INFINITY;
    double sum = 0.0;

    double *TMP = matrix_transpose(A, rows, columns);
    
    for (int i = 0; i < columns; i++) {
	for (int j = 0; j < rows; j++) {
	    sum += TMP[i * rows + j];
	}
	if (max < sum) max = sum;
	sum = 0.0;
    }

    return max;
    
}

double frobenius_norm(double *A, int rows, int columns) {

    double frobenius = 0.0;

    for (int i = 0; i < rows; i++)
	for (int j = 0; j < columns; j++)
	    frobenius += pow(A[i * columns + j], 2.0);
    
    return sqrt(frobenius);
    
}

double matrix_determinant(double *A, int rows, int columns) {

    LU *LU_matrices = LU_decomposition(A, rows, columns);

    double determinant = 1.0;

    for (int i = 0; i < rows; i++)
	for (int j = 0; j < columns; j++)
	    if (i == j) determinant *= LU_matrices->L[i * columns + i];

    for (int i = 0; i < rows; i++)
	for (int j = 0; j < columns; j++)
	    if (i == j) determinant *= LU_matrices->U[i * columns + i];

    free(LU_matrices->L);
    free(LU_matrices->U);
    free(LU_matrices);
    
    return determinant;

}

double *matrix_eigenvalues(double *A, int rows, int columns) {

    // A (rows * columns) ||| Q (rows * rows) ||| R (rows * columns)
    // A (m * n) ||| Q (m * m) ||| R (m * n)
    // m > n with rank(A) = n    
    
    QR *QR_matrices = QR_decomposition(A, rows, columns);

    double min = -INFINITY;

    if (rows < columns)
	min = rows;
    else
	min = columns;

    int index = 0;
    double *eigenvalues = malloc(min * sizeof(double));

    for (int i = 0; i < rows; i++) {
	for (int j = 0; j < columns; j++) {
	    if (i == j) {
		eigenvalues[index] = QR_matrices->R[i * columns + j];
		index++;
	    }
	}
    }

    free(QR_matrices->Q);
    free(QR_matrices->R);
    free(QR_matrices);
    
    return eigenvalues;
    
}
