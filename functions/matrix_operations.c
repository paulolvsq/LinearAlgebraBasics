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

double *forward_substitution(double *L, int d, double *b) {

    double *c = malloc(d * sizeof(double));
    for (int i = 0; i < d; i++)
	c[i] = 0.0;

    for (int i = 0; i < d; i++) {
	c[i] = b[i];
	for (int j = 0; j < i - 1; j++) {
	    c[i] = c[i] - L[i * d + j] * c[j];
	}
	c[i] = c[i] / L[i * d + i];
    }

    return c;

}

double *backward_substitution(double *U, int d, double *c) {

    double *x = malloc(d * sizeof(double));
    for (int i = 0; i < d; i++)
	x[i] = 0.0;

    for (int i = d; i > 0; i--) {
	x[i] = c[i];
	for (int j = i + 1; j < d; j++) {
	    x[i] = x[i] - U[i * d + j] * x[j];
	}
	x[i] = x[i] / U[i * d + i];
    }

    return x;
    
}

double *matrix_inverse(double *A, int d) {

    double *I = generate_identity_matrix(n);

    LU *LU_matrix = LU_decomposition(A, d, d);

    double *C = malloc(d * d * sizeof(double));
    double *inverse = malloc(d * d * sizeof(double));

    double *b_i = malloc(d * sizeof(double));
    
    for (int i = 0; i < d; i++) {
	for (int j = 0; j < d; j++) {
	    b_i[j * d] = I[j * d];
	
    
}

//double *matrix_inverse(double *A, int n) {

    // n by n square matrix A is called invertible if there exists an n by n square matrix B
    // such that AB = BA = I_n where I_n denotes the n by n identity matrix
    // inversible => square matrix and non singular matrix, must have a determinant =/= 0
    /*
    function forward_substitution(matrix L, vector b):
    dimension = size of L
    Initialize vector c of size dimension to zeros

    for i from 1 to dimension:
        c[i] = b[i]
        for j from 1 to i - 1:
            c[i] = c[i] - L[i][j] * c[j]
        c[i] = c[i] / L[i][i]

    return c

function backward_substitution(matrix U, vector c):
    dimension = size of U
    Initialize vector x of size dimension to zeros

    for i from dimension downto 1:
        x[i] = c[i]
        for j from i + 1 to dimension:
            x[i] = x[i] - U[i][j] * x[j]
        x[i] = x[i] / U[i][i]

    return x

function inverse_matrix_with_LU_decomposition(matrix L, matrix U):
    dimension = size of L or U
    Create an identity matrix I of size dimension x dimension

    // Step 1: Résoudre Lc = b pour chaque colonne de la matrice identité
    for i from 1 to dimension:
        Create vector b_i representing the i-th column of the identity matrix I
        Solve the lower triangular system Lc = b_i using forward substitution to find c_i
        Store c_i as the i-th column of a temporary matrix C

    // Step 2: Résoudre Ux = c pour chaque colonne de la matrice C
    for i from 1 to dimension:
        Create vector c_i representing the i-th column of the matrix C
        Solve the upper triangular system Ux = c_i using backward substitution to find x_i
        Store x_i as the i-th column of the inverse matrix

    return the inverse matrix

    */

//}
