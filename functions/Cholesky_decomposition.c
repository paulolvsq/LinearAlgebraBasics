#include "LinearAlgebraBasics.h"
#include <string.h>

/**
 * @brief Allocates and initializes a Cholesky decomposition structure.
 *
 * This function creates a Cholesky decomposition structure for a given matrix A
 * and initializes its components (A, L, and Lᵀ). The input matrix A is copied into
 * the structure, and the matrices L and Lᵀ are initialized to zero.
 *
 * @param A Pointer to the input square matrix (size: size x size).
 * @param size Dimension of the square matrix A (must be positive).
 *
 * @return Pointer to the Cholesky structure on success, or NULL on failure due to
 *         invalid dimensions or memory allocation errors.
 */

Cholesky *create_Cholesky(double *A, int size) {
    
    Cholesky *Cholesky_decomposition = malloc(sizeof(Cholesky));

    if (!Cholesky_decomposition) {
        fprintf(stderr, "Memory allocation failed for Cholesky structure.\n");
        return NULL;
    }

    Cholesky_decomposition->size = size;

    Cholesky_decomposition->A = malloc(size * size * sizeof(double));
    Cholesky_decomposition->L = calloc(size * size, sizeof(double));
    Cholesky_decomposition->L_t = calloc(size * size, sizeof(double));

    if (!Cholesky_decomposition->A || !Cholesky_decomposition->L || !Cholesky_decomposition->L_t) {
        fprintf(stderr, "Memory allocation failed for matrices in Cholesky structure.\n");
        free_Cholesky(Cholesky_decomposition);
        return NULL;
    }

    memcpy(Cholesky_decomposition->A, A, size * size * sizeof(double));

    return Cholesky_decomposition;

}

/**
 * @brief Performs Cholesky decomposition on a symmetric positive-definite matrix.
 *
 * This function decomposes a symmetric positive-definite matrix A into a lower triangular matrix L
 * such that \( A = L \cdot L^T \). The resulting matrices L and Lᵀ are stored in the Cholesky structure.
 *
 * The function checks if the input matrix is symmetric and positive definite before performing the decomposition.
 *
 * @param A Pointer to the input square matrix (size: size x size).
 * @param size Dimension of the square matrix A (must be positive).
 *
 * @return Pointer to the Cholesky structure containing matrices A, L, and Lᵀ on success,
 *         or NULL on failure due to invalid dimensions, non-symmetric matrix, non-positive definite matrix,
 *         or memory allocation errors.
 */

Cholesky *Cholesky_decomposition(double *A, int size) {

    for (int i = 0; i < size; i++) {
        for (int j = i + 1; j < size; j++) {
            if (A[i * size + j] != A[j * size + i]) {
                fprintf(stderr, "Matrix is not symmetric.\n");
                return NULL;
            }
        }
        if (A[i * size + i] <= 0) {
            fprintf(stderr, "Matrix is not positive definite.\n");
            return NULL;
        }
    }

    Cholesky *Cholesky_decomp = create_Cholesky(A, size);
    
    if (!Cholesky_decomp) return NULL;

    for (int i = 0; i < size; i++) {
        for (int j = 0; j <= i; j++) {
            double sum = 0.0;

            if (j == i) { 
                for (int k = 0; k < j; k++)
                    sum += Cholesky_decomp->L[j * size + k] * Cholesky_decomp->L[j * size + k];

                double diag_value = A[j * size + j] - sum;
		
                if (diag_value <= 0) {
                    fprintf(stderr, "Matrix is not positive definite at row %d.\n", j);
                    free_Cholesky(Cholesky_decomp);
                    return NULL;
                }

                Cholesky_decomp->L[j * size + j] = sqrt(diag_value);
            }
	    else { 
                for (int k = 0; k < j; k++)
                    sum += Cholesky_decomp->L[i * size + k] * Cholesky_decomp->L[j * size + k];

                Cholesky_decomp->L[i * size + j] = (A[i * size + j] - sum) / Cholesky_decomp->L[j * size + j];
            }
        }
    }

    double *L_t = matrix_transpose(Cholesky_decomp->L, size, size);
    
    if (!L_t) {
        fprintf(stderr, "Memory allocation failed for transpose of L.\n");
        free_Cholesky(Cholesky_decomp);
        return NULL;
    }

    memcpy(Cholesky_decomp->L_t, L_t, size * size * sizeof(double));
    
    free(L_t);

    return Cholesky_decomp;
}

/**
 * @brief Frees all memory associated with a Cholesky decomposition structure.
 *
 * This function releases the memory allocated for the matrices A, L, and Lᵀ,
 * as well as the Cholesky structure itself.
 *
 * @param Cholesky_decomposition Pointer to the Cholesky structure to free.
 */

void free_Cholesky(Cholesky *Cholesky_decomposition) {

    if (!Cholesky_decomposition) return;

    free(Cholesky_decomposition->A);
    free(Cholesky_decomposition->L);
    free(Cholesky_decomposition->L_t);

    free(Cholesky_decomposition);

}
    
