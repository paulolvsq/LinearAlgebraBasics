#include "LinearAlgebraBasics.h"

double *vectors_addition(double *X, double *Y, int dimension) {

    if (dimension <= 0) {
        fprintf(stderr, "Error: Invalid dimension (%d). Must be strictly positive.\n", dimension);
        return NULL;
    }

    if (!X || !Y) {
        fprintf(stderr, "Error: Null pointer detected in vectors_addition.\n");
        return NULL;
    }

    double *vector = malloc(dimension * sizeof(double));

    if (!vector) {
        fprintf(stderr, "Error: Memory allocation failed for vectors_addition (dimension=%d).\n", dimension);
        return NULL;
    }

    for (int i = 0; i < dimension; i++)
	vector[i] = X[i] + Y[i];

    return vector;
    
}

double *scalar_product(double *X, double *Y, int dimension) {

    if (dimension <= 0) {
        fprintf(stderr, "Error: Invalid dimension (%d). Must be strictly positive.\n", dimension);
        return NULL;
    }

    if (!X || !Y) {
        fprintf(stderr, "Error: Null pointer detected in scalar_product.\n");
        return NULL;
    }

    double *vector = malloc(dimension * sizeof(dimension));

    if (!vector) {
        fprintf(stderr, "Error: Memory allocation failed for scalar_product (dimension=%d).\n", dimension);
        return NULL;
    }
    
    for (int i = 0; i < dimension; i++)
	vector[i] = X[i] * Y[i];

    return vector;

}

double *vector_product(double *X, double *Y) {

    // ONLY IN DIMENSION 3

    const int dimension = 3;

    if (!X || !Y) {
        fprintf(stderr, "Error: Null pointer detected in vector_product.\n");
        return NULL;
    }

    double *vector = malloc(dimension * sizeof(double));

    if (!vector) {
        fprintf(stderr, "Error: Memory allocation failed for vector_product.\n");
        return NULL;
    }

    vector[0] = X[1] * Y[2] - X[2] * Y[1];
    vector[1] = X[2] * Y[0] - X[0] * Y[2];
    vector[2] = X[0] * Y[1] - X[1] * Y[0];

    return vector;
       
}

double vector_norm(double *X, int dimension) {

    if (dimension <= 0) {
        fprintf(stderr, "Error: Invalid dimension (%d). Must be strictly positive.\n", dimension);
	return -1.0; 
    }

    if (!X) {
        fprintf(stderr, "Error: Null pointer detected in vector_norm.\n");
        return -1.0; 
    }

    double sum = 0.0;

    for (int i = 0; i < dimension; i++)
	sum += X[i] * X[i];

    return sqrt(sum);

}
