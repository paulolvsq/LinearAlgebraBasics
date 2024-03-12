#ifndef __LinearAlgebraBasics_
#define __LinearAlgebraBasics_

#include <stdio.h>
#include <stdlib.h>
#include <time.h>

#include <math.h>

/* generate_matrix.c */
double *generate_matrix_double (int rows, int columns);
double *generate_identity_matrix (int dimension);

/* sequential_matrix_product.c */
double *sequential_matrix_product(double *P, int P_rows, int P_columns, double *Q, int Q_rows, int Q_columns);

/* sequential_vector_matrix_product.c */
double *sequential_vector_matrix_product(double *A, int A_rows, int A_columns, double *X, double dimension); 

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

#endif 
