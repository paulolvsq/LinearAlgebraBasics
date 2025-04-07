#include "LinearAlgebraBasics.h"
#include <string.h>

/**
 * @brief Allocates and initializes an LDLT decomposition structure.
 *
 * This function creates an LDLT decomposition structure for a given matrix A
 * and initializes its components (A, L, D, and Lᵀ). The input matrix A is copied
 * into the structure, and the matrices L, D, and Lᵀ are initialized to zero.
 *
 * @param A Pointer to the input square matrix (size: size x size).
 * @param size Dimension of the square matrix A (must be positive).
 *
 * @return Pointer to the LDLT structure on success, or NULL on failure due to
 *         invalid dimensions or memory allocation errors.
 */

LDLT *create_LDLT(double *A, int size) {

    if (!A) {
	fprintf(stderr, "Error: Null pointer detected for input matrix A in create_LDLT.\n");
	return NULL;
    }

    if (size <= 0) {
	fprintf(stderr, "Error: Invalid size (%d). Must be strictly positive.\n", size);
	return NULL;
    }

    LDLT *LDLT_decomposition = malloc(sizeof(LDLT));

    if (!LDLT_decomposition) {
	fprintf(stderr, "Memory allocation failed for LDLT structure.\n");
	return NULL;
    }

    LDLT_decomposition->size = size;

    LDLT_decomposition->A = malloc(size * size * sizeof(double));
    LDLT_decomposition->L = calloc(size * size, sizeof(double));
    LDLT_decomposition->D = calloc(size * size, sizeof(double));
    LDLT_decomposition->L_t = calloc(size * size, sizeof(double));

    if (!LDLT_decomposition->A || !LDLT_decomposition->L || !LDLT_decomposition->D || !LDLT_decomposition->L_t) {
        fprintf(stderr, "Memory allocation failed for matrices in LDLT structure.\n");
        free_LDLT(LDLT_decomposition);
        return NULL;
    }

    memcpy(LDLT_decomposition->A, A, size * size * sizeof(double));

    return LDLT_decomposition;

}

/**
 * @brief Performs LDLT decomposition on a symmetric square matrix.
 *
 * This function decomposes a symmetric square matrix A into:
 * - L: a lower triangular matrix with unit diagonal,
 * - D: a diagonal matrix,
 * - Lᵀ: the transpose of the lower triangular matrix.
 *
 * The decomposition satisfies \( A = L \cdot D \cdot L^T \).
 *
 * The function checks if the input matrix is symmetric before performing the decomposition.
 *
 * @param A Pointer to the input square matrix (size: size x size).
 * @param size Dimension of the square matrix A (must be positive).
 *
 * @return Pointer to the LDLT structure containing matrices A, L, D, and Lᵀ on success,
 *         or NULL on failure due to invalid dimensions, non-symmetric matrix,
 *         singularity detection, or memory allocation errors.
 */

LDLT *LDLT_decomposition(double *A, int size) {

    if (!A) {
        fprintf(stderr, "Error: Null pointer detected for input matrix A in LDLT_decomposition.\n");
        return NULL;
    }

    if (size <= 0) {
        fprintf(stderr, "Error: Invalid size (%d). Must be strictly positive.\n", size);
        return NULL;
    }

    for (int i = 0; i < size; i++) {
	for (int j = i + 1; j < size; j++) {
	    if (A[i * size + j] != A[j * size + i]) {
		fprintf(stderr, "Matrix is not symmetric.\n");
		return NULL;
	    }
	}
    }

    LDLT *LDLT_decomp = create_LDLT(A, size);

    if (!LDLT_decomp) {
	fprintf(stderr, "Memory allocation failed for LDLT structure.\n");
	return NULL;
    }

    for (int i = 0; i < size; i++) {
	LDLT_decomp->L[i * size + i] = 1.0;
    }

    double epsilon = 1e-12;
    
    for (int i = 0; i < size; i++) {
	
	double sum = 0.0;

	for (int k = 0; k < i; k++) {
	    sum += LDLT_decomp->L[i * size + k] * LDLT_decomp->L[i * size + k] * LDLT_decomp->D[k * size + k];
	}

	LDLT_decomp->D[i * size + i] = LDLT_decomp->A[i * size + i] - sum;

	if (fabs(LDLT_decomp->D[i * size + i]) < epsilon) {
	    fprintf(stderr, "Matrix is nearly singular.\n");
	    free_LDLT(LDLT_decomp);
	    return NULL;
	}

	for (int j = i + 1; j < size; j++) {

	    double sum = 0.0;

	    for (int k = 0; k < i; k++) {
		sum += LDLT_decomp->L[j * size + k] * LDLT_decomp->L[i * size + k] * LDLT_decomp->D[k * size + k];
	    }

	    LDLT_decomp->L[j * size + i] = (LDLT_decomp->A[j * size + i] - sum) / LDLT_decomp->D[i * size + i];

	}

    }

    for (int i = 0; i < size; i++) {
	for (int j = 0; j < size; j++) {
	    LDLT_decomp->L_t[j * size + i] = LDLT_decomp->L[i * size + j];
	}
    }
    
    return LDLT_decomp;
    
}

/**
 * @brief Frees all memory associated with an LDLT decomposition structure.
 *
 * This function releases the memory allocated for the matrices A, L, D, and Lᵀ,
 * as well as the LDLT structure itself.
 *
 * @param LDLT_decomposition Pointer to the LDLT structure to free.
 */

void free_LDLT(LDLT *LDLT_decomposition) {

    if (!LDLT_decomposition) return;

    free(LDLT_decomposition->A);
    free(LDLT_decomposition->L);
    free(LDLT_decomposition->L_t);
    free(LDLT_decomposition->D);
    
    free(LDLT_decomposition);

}

