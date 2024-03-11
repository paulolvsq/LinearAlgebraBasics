#include "LinearAlgebraBasics.h"
#define DOUBLE 10.0

double *generate_matrix_double (int rows, int columns) {

    double *matrix = malloc((rows * columns) * sizeof(double)); 
    
    srand(time(NULL));

    for (int i = 0; i < rows; i++) 
	for (int j = 0; j < columns; j++) 
	    matrix[i * columns + j] = (double) (rand()) / RAND_MAX * DOUBLE;    

     return matrix;

}

double *generate_identity_matrix (int dimension) {

    double *matrix = malloc((dimension * dimension) * sizeof(double));

    for (int i = 0; i < dimension; i++) {
	for (int j = 0; j < dimension; j++) {
	    if (i == j) matrix[i * dimension + j] = 1.0;
	    else matrix[i * dimension + j] = 0.0;
	}
    }
    
    return matrix;
    
}
	    
