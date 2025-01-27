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

    free(TMP);
    
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

int has_converged(double *H, int n, double tol) {

    // Suppose that we only use square matrix

    for (int i = 1; i < n; i++) {
        for (int j = 0; j < i; j++) {
            if (fabs(H[i * n + j]) > tol) {
                return 0;
            }
        }
    }
    
    return 1;
		 
}

double *matrix_eigenvalues(double *A, int rows, int columns, int max_iter, double tol) {

    // A (rows * columns) ||| Q (rows * rows) ||| R (rows * columns)
    // A (m * n) ||| Q (m * m) ||| R (m * n)
    // m > n with rank(A) = n    
    
    double min = -INFINITY;

    if (rows < columns)
	min = rows;
    else
	min = columns;

    int size = (int) (min);
    
    double *H = malloc(rows * columns * sizeof(double));
    
    for (int i = 0; i < rows * columns; i++) {
	H[i] = A[i];
    }

    int iter = 0;

    QR *QR_H;
    
    while (iter < max_iter && !has_converged(H, rows, tol)) {

	QR_H = QR_decomposition(H, rows, columns);

	for (int i = 0; i < rows; i++) {
            for (int j = 0; j < columns; j++) {
                H[i * columns + j] = 0.0;
                for (int k = 0; k < columns; k++) {
                    H[i * columns + j] += QR_H->R[i * columns + k] * QR_H->Q[k * columns + j];
                }
            }
        }
        iter++;
	
    }

    double *eigenvalues = malloc(size * sizeof(double));
    
    for (int i = 0; i < size; i++) {
        eigenvalues[i] = H[i * size + i];
    }
    
    printf("Converged in %d iteration(s).\n", iter);

    free(H);
        
    return eigenvalues;
    
}

double *forward_substitution(double *L, int d, double *b) {

    double *c = malloc(d * sizeof(double));

    for (int i = 0; i < d; i++)
	c[i] = 0.0;

    for (int i = 0; i < d; i++) {
	c[i] = b[i];
	for (int j = 0; j < i; j++) {
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

    for (int i = d - 1; i >= 0; i--) {
	x[i] = c[i];
	for (int j = i + 1; j < d; j++) {
	    x[i] = x[i] - U[i * d + j] * x[j];
	}
	x[i] = x[i] / U[i * d + i];
    }

    return x;
    
}

double *solve_LU_system(double *A, double *b, int d) {

    LU *LU_matrix = LU_decomposition(A, d, d);

    // Ly = b
    double *y = forward_substitution(LU_matrix->L, d, b);
    // Ux = y
    double *x = backward_substitution(LU_matrix->U, d, y);

    LU_free(LU_matrix);
    
    free(y);
    
    return x;
    
}

double *matrix_inverse(double *A, int d) {

    double *I = generate_identity_matrix(d);

    LU *LU_matrix = LU_decomposition(A, d, d);

    double *C = malloc(d * d * sizeof(double));
    double *inverse = malloc(d * d * sizeof(double));

    double *b_i = malloc(d * sizeof(double));

    double *c_i;
    double *x_i;
    
    for (int i = 0; i < d; i++) {
	for (int j = 0; j < d; j++) {
	    b_i[j] = I[j * d + i];
	}
	c_i = forward_substitution(LU_matrix->L, d, b_i);
	for (int j = 0; j < d; j++) {
	    C[j * d + i] = c_i[j];
	}
	free(c_i);
    }

    for (int i = 0; i < d; i++) {
	for (int j = 0; j < d; j++) {
	    b_i[j] = C[j * d + i];
	}
	x_i = backward_substitution(LU_matrix->U, d, b_i);
	for (int j = 0; j < d; j++) {
	    inverse[j * d + i] = x_i[j];
	}
	free(x_i);
    }

    free(I);
    free(C);
    free(b_i);

    LU_free(LU_matrix);
    
    return inverse;
    
}
