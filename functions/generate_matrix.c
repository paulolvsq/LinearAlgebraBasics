#include "LinearAlgebraBasics.h"
#define INTEGER 10
#define DOUBLE 10.0

double *generate_matrix_double (int rows, int columns) {

    double *matrix = malloc((rows * columns) * sizeof(double)); 
    
    srand(time(NULL));

    for (int i = 0; i < rows; i++) 
	for (int j = 0; j < columns; j++) 
	    matrix[i * columns +j] = (double) (rand()) / RAND_MAX * DOUBLE;    

     return matrix;

}


