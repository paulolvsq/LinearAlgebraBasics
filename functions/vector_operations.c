#include "LinearAlgebraBasics.h"

/**
 * @brief Computes the element-wise addition of two vectors X and Y.
 *
 * This function calculates the vector result = X + Y, where X and Y are vectors
 * of the same dimension.
 *
 * @param X Pointer to the first input vector (size: dimension).
 * @param Y Pointer to the second input vector (size: dimension).
 * @param dimension Dimension of the vectors (must be positive).
 *
 * @return Pointer to the resulting vector on success, or NULL on failure due to invalid dimension,
 *         null pointers, or memory allocation errors.
 */

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

/**
 * @brief Computes the element-wise product of two vectors X and Y.
 *
 * This function calculates the vector result = X .* Y, where .* denotes
 * element-wise multiplication.
 *
 * @param X Pointer to the first input vector (size: dimension).
 * @param Y Pointer to the second input vector (size: dimension).
 * @param dimension Dimension of the vectors (must be positive).
 *
 * @return Pointer to the resulting vector on success, or NULL on failure due to invalid dimension,
 *         null pointers, or memory allocation errors.
 */

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

/**
 * @brief Computes the cross product of two 3-dimensional vectors X and Y.
 *
 * This function calculates the vector result = X x Y, where x denotes
 * the cross product. It assumes that both input vectors are 3-dimensional.
 *
 * @param X Pointer to the first input vector (size: 3).
 * @param Y Pointer to the second input vector (size: 3).
 *
 * @return Pointer to the resulting vector (size: 3) on success, or NULL on failure due to null pointers
 *         or memory allocation errors.
 */

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

/**
 * @brief Computes the Euclidean norm (length) of a given vector X.
 *
 * This function calculates the norm ||X|| = sqrt(X_1^2 + X_2^2 + ... + X_n^2),
 * where n is the dimension of the vector.
 *
 * @param X Pointer to the input vector (size: dimension).
 * @param dimension Dimension of the vector (must be positive).
 *
 * @return The Euclidean norm of the vector on success, or -1.0 on failure due to invalid dimension
 *         or null pointers.
 */

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
