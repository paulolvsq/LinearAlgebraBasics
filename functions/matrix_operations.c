#include "LinearAlgebraBasics.h"

// ASSERTION : THE USER WILL ALWAYS GIVE MATRICES OF THE RIGHT SIZE AS INPUT 

double *matrices_addition(double *A, double *B, int rows, int columns) {

    if (rows <= 0 || columns <= 0) {
        fprintf(stderr, "Error: Invalid dimensions for matrix addition (rows=%d, columns=%d). Both must be strictly positive.\n", rows, columns);
        return NULL;
    }

    if (!A || !B) {
        fprintf(stderr, "Error: Null pointer detected for input matrices in matrices_addition.\n");
        return NULL;
    }

    double *matrix = malloc(rows * columns * sizeof(double));

    if (!matrix) {
        fprintf(stderr, "Error: Memory allocation failed for result matrix in matrices_addition.\n");
        return NULL;
    }

    for (int i = 0; i < rows; i++)
	for (int j = 0; j < columns; j++)
	    matrix[i * columns + j] = A[i * columns + j] + B[i * columns + j];

    return matrix;

}

double *matrix_scalar_multiplication(double *A, int rows, int columns, int scalar) {

    if (rows <= 0 || columns <= 0) {
        fprintf(stderr, "Error: Invalid dimensions for scalar multiplication (rows=%d, columns=%d). Both must be strictly positive.\n", rows, columns);
        return NULL;
    }

    double *matrix = malloc(rows * columns * sizeof(double));

    if (!A) {
        fprintf(stderr, "Error: Null pointer detected for input matrix in matrix_scalar_multiplication.\n");
        return NULL;
    }

    if (!matrix) {
        fprintf(stderr, "Error: Memory allocation failed for result matrix in matrix_scalar_multiplication.\n");
        return NULL;
    }
    
    for (int i = 0; i < rows; i++)
	for (int j = 0; j < columns; j++)
	    matrix[i * columns + j] = scalar * A[i * columns + j];

    return matrix;
	       
}

double *matrix_transpose(double *A, int rows, int columns) {

    if (rows <= 0 || columns <= 0) {
        fprintf(stderr, "Error: Invalid dimensions for transpose (rows=%d, columns=%d). Both must be strictly positive.\n", rows, columns);
        return NULL;
    }

    if (!A) {
        fprintf(stderr, "Error: Null pointer detected for input matrix in matrix_transpose.\n");
        return NULL;
    }

    double *matrix = malloc(rows * columns * sizeof(double));

    if (!matrix) {
        fprintf(stderr, "Error: Memory allocation failed for transposed matrix.\n");
        return NULL;
    }
    
    for (int i = 0; i < columns; i++)
	for (int j = 0; j < rows; j++)
	    matrix[i * rows + j] = A[j * rows + i];

    return matrix;

}

double matrix_trace(double *A, int dimension) {

    if (dimension <= 0) {
        fprintf(stderr, "Error: Invalid dimension (%d). Must be strictly positive.\n", dimension);
        return -1.0;
    }

    if (!A) {
        fprintf(stderr, "Error: Null pointer detected for input square matrix in matrix_trace.\n");
        return -1.0;
    }

    double trace = 0.0;

    for (int i = 0; i < dimension; i++) 
        trace += A[i * dimension + i];
	
    return trace;
    
}

double matrix_norm(double *A, int rows, int columns) {

    // 1-NORM

    if (rows <= 0 || columns <= 0) {
        fprintf(stderr, "Error: Invalid dimensions for matrix norm computation (rows=%d, columns=%d). Both must be strictly positive.\n", rows, columns);
        return -1.0; 
    }

    if (!A) {
        fprintf(stderr, "Error: Null pointer detected for input matrix in matrix_norm.\n");
        return -1.0;
    }

    double max = -INFINITY;
    double sum = 0.0;

    double *TMP = matrix_transpose(A, rows, columns);

    if (!TMP) {
        fprintf(stderr, "Error: Failed to transpose matrix in matrix_norm.\n");
        return -1.0;
    }
    
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

    if (rows <= 0 || columns <= 0) {
        fprintf(stderr, "Error: Invalid dimensions for Frobenius norm computation (rows=%d, columns=%d). Both must be strictly positive.\n", rows, columns);
        return -1.0;
    }

    if (!A) {
        fprintf(stderr, "Error: Null pointer detected for input matrix in Frobenius norm computation.\n");
        return -1.0; 
    }

    double frobenius = 0.0;

    for (int i = 0; i < rows; i++)
	for (int j = 0; j < columns; j++)
	    frobenius += pow(A[i * columns + j], 2.0);
    
    return sqrt(frobenius);
    
}

double matrix_determinant(double *A, int rows, int columns) {

    if (rows <= 0 || columns <= 0 || rows != columns) {
        fprintf(stderr, "Error: Invalid dimensions for determinant computation (rows=%d, columns=%d). Must be a square matrix.\n", rows, columns);
        return -1.0; 
    }

    if (!A) {
        fprintf(stderr, "Error: Null pointer detected for input square matrix in determinant computation.\n");
        return -1.0; 
    }
    
    LU *LU_matrices = LU_decomposition(A, rows, columns);

    if (!LU_matrices) {
        fprintf(stderr, "Error: LU decomposition failed during determinant computation.\n");
        return -1.0; 
    }

    double determinant = 1.0;
    double epsilon = 1e-10;
    
    for (int i = 0; i < rows; i++) {
	
        determinant *= LU_matrices->L[i * columns + i];
        determinant *= LU_matrices->U[i * columns + i];
	
        if (fabs(LU_matrices->L[i * columns + i]) < epsilon || fabs(LU_matrices->U[i * columns + i]) < epsilon) {
            fprintf(stderr, "Error: Singular matrix detected during determinant computation.\n");
            LU_free(LU_matrices);
            return -1.0;
        }
	
    }
    
    LU_free(LU_matrices);
    
    return determinant;

}

int has_converged(double *H, int n, double tol) {

    if (n <= 0) {
        fprintf(stderr, "Error: Invalid dimension (%d). Must be strictly positive.\n", n);
        return -1; 
    }

    if (tol <= 0) {
        fprintf(stderr, "Error: Invalid tolerance (%f). Must be strictly positive.\n", tol);
        return -1; 
    }

    if (!H) {
        fprintf(stderr, "Error: Null pointer detected for input matrix in has_converged.\n");
        return -1; 
    }

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

    if (rows <= 0 || columns <= 0) {
        fprintf(stderr, "Error: Invalid dimensions for eigenvalue computation (rows=%d, columns=%d). Both must be strictly positive.\n", rows, columns);
        return NULL;
    }

    if (tol <= 0) {
        fprintf(stderr, "Error: Invalid tolerance (%f). Must be strictly positive.\n", tol);
        return NULL;
    }

    if (!A) {
        fprintf(stderr, "Error: Null pointer detected for input matrix in matrix_eigenvalues.\n");
        return NULL;
    }
    
    double min = rows < columns ? rows : columns;
    int size = (int) (min);
    
    double *H = malloc(rows * columns * sizeof(double));

    if (!H) {
        fprintf(stderr, "Error: Memory allocation failed for intermediate matrix H.\n");
        return NULL;
    }
    
    for (int i = 0; i < rows * columns; i++) {
	H[i] = A[i];
    }

    int iter = 0;

    QR *QR_H;
    
    while (iter < max_iter && !has_converged(H, rows, tol)) {

	QR_H = QR_decomposition(H, rows, columns);

	if (!QR_H) {
            fprintf(stderr, "Error: QR decomposition failed during eigenvalue computation.\n");
            free(H);
            return NULL;
        }

	for (int i = 0; i < rows; i++) {
            for (int j = 0; j < columns; j++) {
                H[i * columns + j] = 0.0;
                for (int k = 0; k < columns; k++) {
                    H[i * columns + j] += QR_H->R[i * columns + k] * QR_H->Q[k * columns + j];
                }
            }
        }
        iter++;

	if (iter >= max_iter) {
	    fprintf(stderr, "Maximum iteration reached without convergence\n");
	    free(H);
	    return NULL;
	}
	
    }

    QR_free(QR_H);

    
    double *eigenvalues = malloc(size * sizeof(double));

    if (!eigenvalues) {
	fprintf(stderr,"Allocation failed\n");
	free(H);
	return NULL;
    }
    
    for (int i = 0; i < size; i++) {
        eigenvalues[i] = H[i * size + i];
    }
    
    printf("Converged in %d iteration(s).\n", iter);

    free(H);
        
    return eigenvalues;
    
}

double *forward_substitution(double *L, int d, double *b) {

    if (d <= 0) {
        fprintf(stderr, "Error: Invalid dimension (%d). Must be strictly positive.\n", d);
        return NULL;
    }

    if (!L || !b) {
        fprintf(stderr, "Error: Null pointer detected in forward_substitution.\n");
        return NULL;
    }

    double *c = malloc(d * sizeof(double));

    if (!c) {
        fprintf(stderr, "Error: Memory allocation failed for solution vector in forward_substitution.\n");
        return NULL;
    }

    double epsilon = 1e-10;

    for (int i = 0; i < d; i++)
	c[i] = 0.0;

    for (int i = 0; i < d; i++) {
	if (fabs(L[i * d + i]) < epsilon) { 
            fprintf(stderr, "Error: Singular matrix detected in forward_substitution at row %d.\n", i);
            free(c);
            return NULL;
        }
	c[i] = b[i];
	for (int j = 0; j < i; j++) {
	    c[i] = c[i] - L[i * d + j] * c[j];
	}
	c[i] = c[i] / L[i * d + i];
    }

    return c;

}

double *backward_substitution(double *U, int d, double *c) {

    if (d <= 0) {
        fprintf(stderr, "Error: Invalid dimension (%d). Must be strictly positive.\n", d);
        return NULL;
    }

    if (!U || !c) {
        fprintf(stderr, "Error: Null pointer detected in backward_substitution.\n");
        return NULL;
    }

    double *x = malloc(d * sizeof(double));

    if (!x) {
        fprintf(stderr, "Error: Memory allocation failed for solution vector in backward_substitution.\n");
        return NULL;
    }

    for (int i = 0; i < d; i++)
	x[i] = 0.0;

    for (int i = d - 1; i >= 0; i--) {
	if (fabs(U[i * d + i]) < 1e-10) { 
            fprintf(stderr, "Error: Singular matrix detected in backward_substitution at row %d.\n", i);
            free(x);
            return NULL;
        }
	x[i] = c[i];
	for (int j = i + 1; j < d; j++) {
	    x[i] = x[i] - U[i * d + j] * x[j];
	}
	x[i] = x[i] / U[i * d + i];
    }

    return x;
    
}

double *solve_LU_system(double *A, double *b, int d) {

    if (d <= 0) {
        fprintf(stderr, "Error: Invalid dimension (%d). Must be strictly positive.\n", d);
        return NULL;
    }

    if (!A || !b) {
        fprintf(stderr, "Error: Null pointer detected in solve_LU_system.\n");
        return NULL;
    }

    LU *LU_matrix = LU_decomposition(A, d, d);

    if (!LU_matrix) {
        fprintf(stderr, "Error: LU decomposition failed in solve_LU_system.\n");
        return NULL;
    }
    
    // Ly = b
    double *y = forward_substitution(LU_matrix->L, d, b);

    if (!y) {
        fprintf(stderr, "Error: Forward substitution failed in solve_LU_system.\n");
        LU_free(LU_matrix);
        return NULL;
    }
    
    // Ux = y
    double *x = backward_substitution(LU_matrix->U, d, y);

    if (!x) {
        fprintf(stderr, "Error: Backward substitution failed in solve_LU_system.\n");
        free(y);
        LU_free(LU_matrix);
        return NULL;
    }

    LU_free(LU_matrix);
    
    free(y);
    
    return x;
    
}

double *matrix_inverse(double *A, int d) {

    if (d <= 0) {
        fprintf(stderr, "Error: Invalid dimension (%d). Must be strictly positive.\n", d);
        return NULL;
    }

    if (!A) {
        fprintf(stderr, "Error: Null pointer detected for input matrix in matrix_inverse.\n");
        return NULL;
    }

    double *I = generate_identity_matrix(d);

    if (!I) {
        fprintf(stderr, "Error: Failed to generate identity matrix in matrix_inverse.\n");
        return NULL;
    }

    LU *LU_matrix = LU_decomposition(A, d, d);

    if (!LU_matrix) {
        fprintf(stderr, "Error: LU decomposition failed in matrix_inverse.\n");
        free(I);
        return NULL;
    }

    double *C = malloc(d * d * sizeof(double));
    double *inverse = malloc(d * d * sizeof(double));
    double *b_i = malloc(d * sizeof(double));

    if (!C || !inverse || !b_i) {
        fprintf(stderr, "Error: Memory allocation failed for matrices in matrix_inverse.\n");
        free(I);
        free(C);
        free(inverse);
        free(b_i);
        LU_free(LU_matrix);
        return NULL;
    }

    double *c_i;
    double *x_i;
    
    for (int i = 0; i < d; i++) {

	for (int j = 0; j < d; j++) {
	    b_i[j] = I[j * d + i];
	}
	
	c_i = forward_substitution(LU_matrix->L, d, b_i);

	if (!c_i) {
            fprintf(stderr, "Error: Forward substitution failed in matrix_inverse.\n");
            free(I);
            free(C);
            free(inverse);
            free(b_i);
            LU_free(LU_matrix);
            return NULL;

	}
	
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

	if (!x_i) {
	    fprintf(stderr,"Backward substitution failed\n");
	    free(I);
	    free(C);
	    free(b_i);
	    LU_free(LU_matrix);
	    return NULL;

	}
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
