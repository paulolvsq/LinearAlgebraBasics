#ifndef __LinearAlgebraBasics_
#define __LinearAlgebraBasics_

#include <stdio.h>
#include <stdlib.h>
#include <time.h>

#include <math.h>

typedef struct LU {

    int rows, columns;

    double *A;
    double *L;
    double *U;

} LU;

typedef struct QR {

    int rows, columns;

    double *A;
    double *Q;
    double *R;

} QR;

/* generate_matrix.c */
double *generate_matrix_double (int rows, int columns);
double *generate_identity_matrix (int dimension);

/* sequential_matrix_product.c */
double *sequential_matrix_product(double *P, int P_rows, int P_columns, double *Q, int Q_rows, int Q_columns);

/* sequential_vector_matrix_product.c */
double *sequential_vector_matrix_product(double *A, int A_rows, int A_columns, double *X, double dimension); 

/* parallel_vector_matrix_product.c */
double *parallel_vector_matrix_product(double *A, int A_rows, int A_columns, double *X, double dimension);

/* parallel_matrix_product.c */
double *parallel_matrix_product(double *P, int P_rows, int P_columns, double *Q, int Q_rows, int Q_columns);

/* vector_operations.c */
double *vectors_addition(double *X, double *Y, int dimension);
double *scalar_product(double *X, double *Y, int dimension);
double *vector_product(double *X, double *Y);
double vector_norm(double *X, int dimension);

/* matrix_operations.c */
double *matrices_addition(double *A, double *B, int rows, int columns);
double *matrix_scalar_multiplication(double *A, int rows, int columns, int scalar);
double *matrix_transpose(double *A, int rows, int columns);
double matrix_trace(double *A, int dimension);
double matrix_norm(double *A, int rows, int columns);
double frobenius_norm(double *A, int rows, int columns);
double matrix_determinant(double *A, int rows, int columns);
double *matrix_eigenvalues(double *A, int rows, int columns, int max_iter, double tol);
double *forward_substitution(double *L, int d, double *b);
double *backward_substitution(double *U, int d, double *c);
double *solve_LU_system(double *A, double *b, int d);
double *matrix_inverse(double *A, int n);
int has_converged(double *H, int n, double tol);

/* LU_decomposition.c */
LU *create_LU(double *A, int rows, int columns);
LU *LU_decomposition(double *A, int rows, int columns);
LU *LU_decomposition_parallel(double *A, int rows, int columns);
void LU_free(LU *LU_decomposition);

/* QR_decomposition.c */
QR *create_QR(double *A, int rows, int columns);
QR *QR_decomposition(double *A, int rows, int columns);
QR *QR_decomposition_parallel(double *A, int rows, int columns);
void QR_free(QR *QR_decomposition);

#endif 
