#include "LinearAlgebraBasics.h"

LU *create_LU(double *A, int dimension) {

    LU *LU_decomposition = malloc(sizeof(LU));

    LU_decomposition->A = A;
    LU_decomposition->L = malloc(dimension * sizeof(double));
    LU_decomposition->U = malloc(dimension * sizeof(double));
    
    return LU_decomposition;

}

LU *LU_decomposition(double *A, int dimension) {

    // ASSERTION : THE USER WILL ALWAYS GIVE MATRICES OF THE RIGHT SIZE AS INPUT
    
    LU *LU_decomposition = create_LU(A, dimension);

    for (int i = 0; i < dimension; i++) {
	for (int j = 0; j <= i; j++) {
	    if (j < i) {
		LU_decomposition->L[i * dimension + j] = A[i * dimension + j];
	    }
	    else if (j == i) {
		LU_decomposition->L[i * dimension + j] = 1;
	    }
	    else {
		LU_decomposition->L[i * dimension + j] = 0;
	    }
	}
    }

    for (int i = 0; i < dimension; i++) {
	for (int j = i; j < dimension; j++) {
	    if (j < i) {
		LU_decomposition->U[i * dimension + j] = 0;
	    }
	    else {
		LU_decomposition->U[i * dimension + j] = A[i * dimension + j];
	    }
	}
    }

    for (int k = 0; k < dimension - 1; k++) {
	for (int i = k + 1; i < dimension; i++) {
	    LU_decomposition->L[i * dimension + k] = A[i * dimension + k] / LU_decomposition->U[k * dimension + k];
	    for (int j = k; j < dimension; j++) {
		LU_decomposition->U[i * dimension + j] = LU_decomposition->U[i * dimension + j] - LU_decomposition->L[i * dimension + k] * LU_decomposition->U[k * dimension + j] / LU_decomposition->U[k * dimension + k];
	        
	    }
	}
    }

    return LU_decomposition;

}

/* fonction decomposition_LU(A, L, U) */
/*     n <- nombre de lignes de A */
    
/*     pour i de 1 à n faire */
/*         pour j de 1 à n faire */
/*             si j < i alors */
/*                 L[i][j] <- A[i][j] */
/*             sinon si j == i alors */
/*                 L[i][j] <- 1 */
/*             sinon */
/*                 L[i][j] <- 0 */
/*             fin si */
/*         fin pour */
/*     fin pour */
    
/*     pour i de 1 à n faire */
/*         pour j de 1 à n faire */
/*             si j < i alors */
/*                 U[i][j] <- 0 */
/*             sinon */
/*                 U[i][j] <- A[i][j] */
/*             fin si */
/*         fin pour */
/*     fin pour */
    
/*     pour k de 1 à n-1 faire */
/*         pour i de k+1 à n faire */
/*             si U[k][k] == 0 alors */
/*                 afficher "La matrice est singulière, impossible de procéder à la décomposition LU" */
/*                 retourner */
/*             fin si */
/*             L[i][k] <- U[i][k] / U[k][k] */
/*             pour j de k à n faire */
/*                 U[i][j] <- U[i][j] - L[i][k] * U[k][j] */
/*             fin pour */
/*         fin pour */
/*     fin pour */
