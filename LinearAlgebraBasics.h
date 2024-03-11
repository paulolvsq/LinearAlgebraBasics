#ifndef __LinearAlgebraBasics_
#define __LinearAlgebraBasics_

#include <stdio.h>
#include <stdlib.h>
#include <time.h>

/* generate_matrix.c */
double *generate_matrix_double (int rows, int columns);
double *generate_identity_matrix (int dimension);

/* sequential_matrix_product.c */
double *sequential_matrix_product(double *P, int P_rows, int P_columns, double *Q, int Q_rows, int Q_columns);

#endif 
