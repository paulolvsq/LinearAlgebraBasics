#include "LinearAlgebraBasics.h"

LU *create_LU(double *A, int dimension) {

    LU *LU_decomposition = malloc(sizeof(LU));

    LU_decomposition->dimension = dimension;

    LU_decomposition->A = A;
    LU_decomposition->L = malloc(dimension * dimension * sizeof(double));
    LU_decomposition->U = malloc(dimension * dimension * sizeof(double));
    
    return LU_decomposition;

}

LU *LU_decomposition(double *A, int dimension) {

    // ASSERTION : THE USER WILL ALWAYS GIVE MATRICES OF THE RIGHT SIZE AS INPUT
    
    LU *LU_decomposition = create_LU(A, dimension);

    for (int i = 0; i < dimension; i++) {
        for (int j = 0; j < dimension; j++) {
            if (i == j)
                LU_decomposition->L[i * dimension + j] = 1;
            else
                LU_decomposition->L[i * dimension + j] = 0;
            LU_decomposition->U[i * dimension + j] = 0;
        }
    }

    for (int i = 0; i < dimension; i++) {
        for (int j = i; j < dimension; j++) {
            LU_decomposition->U[i * dimension + j] = A[i * dimension + j];
            for (int k = 0; k < i; k++) {
                LU_decomposition->U[i * dimension + j] -= LU_decomposition->L[i * dimension + k] * LU_decomposition->U[k * dimension + j];
            }
        }

        for (int j = i + 1; j < dimension; j++) {
            LU_decomposition->L[j * dimension + i] = A[j * dimension + i];
            for (int k = 0; k < i; k++) {
                LU_decomposition->L[j * dimension + i] -= LU_decomposition->L[j * dimension + k] * LU_decomposition->U[k * dimension + i];
            }
            LU_decomposition->L[j * dimension + i] /= LU_decomposition->U[i * dimension + i];
        }
    }
    
    return LU_decomposition;

}
