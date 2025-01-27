#include "LinearAlgebraBasics.h"
#include <string.h>

LU *create_LU(double *A, int rows, int columns) {

    LU *LU_decomposition = malloc(sizeof(LU));

    LU_decomposition->rows = rows;
    LU_decomposition->columns = columns;

    LU_decomposition->A = malloc(rows * columns * sizeof(double));

    memcpy(LU_decomposition->A, A, rows * columns * sizeof(double));
    
    LU_decomposition->L = malloc(rows * columns * sizeof(double));
    LU_decomposition->U = malloc(rows * columns * sizeof(double));
    
    return LU_decomposition;

}

LU *LU_decomposition(double *A, int rows, int columns) {

    // ASSERTION : THE USER WILL ALWAYS GIVE MATRICES OF THE RIGHT SIZE AS INPUT
    
    LU *LU_decomposition = create_LU(A, rows, columns);

    for (int i = 0; i < rows; i++) {
        for (int j = 0; j < columns; j++) {
            if (i == j)
                LU_decomposition->L[i * columns + j] = 1;
            else
                LU_decomposition->L[i * columns + j] = 0;
            LU_decomposition->U[i * columns + j] = 0;
        }
    }

    for (int i = 0; i < rows; i++) {
        for (int j = i; j < columns; j++) {
            LU_decomposition->U[i * columns + j] = A[i * columns + j];
            for (int k = 0; k < i; k++) {
                LU_decomposition->U[i * columns + j] -= LU_decomposition->L[i * columns + k] * LU_decomposition->U[k * columns + j];
            }
        }

        for (int j = i + 1; j < columns; j++) {
            LU_decomposition->L[j * columns + i] = A[j * columns + i];
            for (int k = 0; k < i; k++) {
                LU_decomposition->L[j * columns + i] -= LU_decomposition->L[j * columns + k] * LU_decomposition->U[k * columns + i];
            }
            LU_decomposition->L[j * columns + i] /= LU_decomposition->U[i * columns + i];
        }
    }
    
    return LU_decomposition;

}

LU *LU_decomposition_parallel(double *A, int rows, int columns) {

    // ASSERTION : THE USER WILL ALWAYS GIVE MATRICES OF THE RIGHT SIZE AS INPUT

        LU *LU_decomposition = create_LU(A, rows, columns);
    
#pragma omp parallel for
    for (int i = 0; i < rows; i++) {
        for (int j = 0; j < columns; j++) {
            if (i == j)
                LU_decomposition->L[i * columns + j] = 1;
            else
                LU_decomposition->L[i * columns + j] = 0;
            LU_decomposition->U[i * columns + j] = 0;
        }
    }

    for (int i = 0; i < rows; i++) {
#pragma omp parallel for
        for (int j = i; j < columns; j++) {
            double sum = 0.0;
            for (int k = 0; k < i; k++) {
                sum += LU_decomposition->L[i * columns + k] * LU_decomposition->U[k * columns + j];
            }
            LU_decomposition->U[i * columns + j] = A[i * columns + j] - sum;
        }
	
#pragma omp parallel for
        for (int j = i + 1; j < rows; j++) {
            double sum = 0.0;
            for (int k = 0; k < i; k++) {
		sum += LU_decomposition->L[j * columns + k] * LU_decomposition->U[k * columns + i];
            }
            LU_decomposition->L[j * columns + i] = (A[j * columns + i] - sum) / LU_decomposition->U[i * columns + i];
        }
    }

    return LU_decomposition;
    
}

void LU_free(LU *LU_decomposition) {

    free(LU_decomposition->A);
    free(LU_decomposition->L);
    free(LU_decomposition->U);

    free(LU_decomposition);

}
