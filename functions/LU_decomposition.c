#include "LinearAlgebraBasics.h"
#include <string.h>

/**
 * @brief Allocates and initializes an LU decomposition structure.
 *
 * This function creates an LU decomposition structure for a given matrix A
 * and initializes its components (L, U, and A).
 *
 * @param A Pointer to the input matrix (size: rows x columns).
 * @param rows Number of rows in the matrix (must be positive).
 * @param columns Number of columns in the matrix (must be positive).
 *
 * @return Pointer to the LU structure containing L, U, and A matrices on success,
 *         or NULL on failure due to invalid dimensions or memory allocation errors.
 */

LU *create_LU(double *A, int rows, int columns) {

    if (rows <= 0 || columns <= 0) {
        fprintf(stderr, "Error: Invalid dimensions for LU decomposition (rows=%d, columns=%d). Both must be strictly positive.\n", rows, columns);
        return NULL;
    }

    if (!A) {
        fprintf(stderr, "Error: Null pointer detected for input matrix in create_LU.\n");
        return NULL;
    }

    LU *LU_decomposition = malloc(sizeof(LU));

    if (!LU_decomposition) {
        fprintf(stderr, "Error: Memory allocation failed for LU structure.\n");
        return NULL;
    }

    LU_decomposition->rows = rows;
    LU_decomposition->columns = columns;

    LU_decomposition->A = malloc(rows * columns * sizeof(double));
    LU_decomposition->L = malloc(rows * columns * sizeof(double));
    LU_decomposition->U = malloc(rows * columns * sizeof(double));

    if (!LU_decomposition->A || !LU_decomposition->L || !LU_decomposition->U) {
        fprintf(stderr, "Error: Memory allocation failed for matrices in create_LU.\n");
        free(LU_decomposition->A);
        free(LU_decomposition->L);
        free(LU_decomposition->U);
        free(LU_decomposition);
        return NULL;
    }

    memcpy(LU_decomposition->A, A, rows * columns * sizeof(double));

    return LU_decomposition;

}

/**
 * @brief Performs LU decomposition on a given matrix using Doolittle's algorithm.
 *
 * This function decomposes a matrix A into a lower triangular matrix L
 * and an upper triangular matrix U such that A = L * U.
 *
 * @param A Pointer to the input matrix (size: rows x columns).
 * @param rows Number of rows in the matrix (must be positive).
 * @param columns Number of columns in the matrix (must be positive).
 *
 * @return Pointer to the LU structure containing L and U matrices on success,
 *         or NULL on failure due to invalid dimensions or memory allocation errors.
 */

LU *LU_decomposition(double *A, int rows, int columns) {
 
    LU *LU_decomposition = create_LU(A, rows, columns);

    if (!LU_decomposition) {
        fprintf(stderr, "Error: Failed to create LU decomposition structure.\n");
        return NULL;
    }

    for (int i = 0; i < rows; i++) {
        for (int j = 0; j < columns; j++) {
            if (i == j) {
                LU_decomposition->L[i * columns + j] = 1;
	    }
            else {
                LU_decomposition->L[i * columns + j] = 0;
	    }
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

/**
 * @brief Performs parallelized LU decomposition on a given matrix using OpenMP.
 *
 * This function decomposes a matrix A into a lower triangular matrix L
 * and an upper triangular matrix U such that A = L * U. The computation is
 * parallelized using OpenMP for improved performance on large matrices.
 *
 * @param A Pointer to the input matrix (size: rows x columns).
 * @param rows Number of rows in the matrix (must be positive).
 * @param columns Number of columns in the matrix (must be positive).
 *
 * @return Pointer to the LU structure containing L and U matrices on success,
 *         or NULL on failure due to invalid dimensions or memory allocation errors.
 */

LU *LU_decomposition_parallel(double *A, int rows, int columns) {

    // ASSERTION : THE USER WILL ALWAYS GIVE MATRICES OF THE RIGHT SIZE AS INPUT

    LU *LU_decomposition = create_LU(A, rows, columns);

    if (!LU_decomposition) {
        fprintf(stderr, "Error: Failed to create LU decomposition structure.\n");
        return NULL;
    }
    
#pragma omp parallel for
    for (int i = 0; i < rows; i++) {
        for (int j = 0; j < columns; j++) {
            if (i == j) {
                LU_decomposition->L[i * columns + j] = 1;
	    }
	    else {
                LU_decomposition->L[i * columns + j] = 0;
	    }
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

/**
 * @brief Frees all memory associated with an LU decomposition structure.
 *
 * This function releases the memory allocated for the matrices L, U, and A,
 * as well as the LU structure itself.
 *
 * @param LU_decomposition Pointer to the LU structure to free.
 */

void LU_free(LU *LU_decomposition) {

    if (!LU_decomposition) return;
    
    free(LU_decomposition->A);
    free(LU_decomposition->L);
    free(LU_decomposition->U);

    free(LU_decomposition);

}
