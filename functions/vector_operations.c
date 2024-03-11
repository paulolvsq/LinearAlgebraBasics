#include "LinearAlgebraBasics.h"
#include <math.h> 

double *vectors_addition(double *X, double *Y, int dimension) {

    double *vector = malloc(dimension * sizeof(double));

    for (int i = 0; i < dimension; i++)
	vector[i] = X[i] + Y[i];

    return vector;
    
}

double *scalar_product(double *X, double *Y, int dimension) {

    double *vector = malloc(dimension * sizeof(dimension));
    
    for (int i = 0; i < dimension; i++)
	vector[i] = X[i] * Y[i];

    return vector;

}

double *vector_product(double *X, double *Y) {

    // ONLY IN DIMENSION 3

    double *vector = malloc(3 * sizeof(double));

    vector[0] = X[1] * Y[2] - X[2] * Y[1];
    vector[1] = X[2] * Y[0] - X[0] * Y[2];
    vector[2] = X[0] * Y[1] - X[1] * Y[0];

    return vector;
       
}

double vector_norm(double *X, int dimension) {

    double sum = 0.0;

    for (int i = 0; i < dimension; i++)
	sum += X[i] * X[i];

    return sqrt(sum);

}
