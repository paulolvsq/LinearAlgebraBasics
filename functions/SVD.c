#include "LinearAlgebraBasics.h"
#include <string.h>

SVD *create_SVD(double *A, int rows, int columns) {

    if (!A) {
	fprintf(stderr, "Error: Null pointer detected for input matrix A in create_SVD.\n");
	return NULL;
    }

    if (rows <= 0 || columns <= 0) {
	fprintf(stderr, "Error: Invalid size in SVD (%d, %d). Must be strictly positive.\n", rows, columns);
	return NULL;
    }

    SVD *SVD_decomposition = malloc(sizeof(SVD));

    if (!SVD_decomposition) {
	fprintf(stderr, "Memory allocation failed for SVD structure.\n");
	return NULL;
    }

    int min = (rows < columns) ? rows : columns;

    SVD_decomposition->rows = rows;
    SVD_decomposition->columns = columns;

    SVD_decomposition->A = malloc(rows * columns * sizeof(double));
    SVD_decomposition->U = calloc(rows * rows, sizeof(double));
    SVD_decomposition->S = calloc(min, sizeof(double));
    SVD_decomposition->V = calloc(columns * columns, sizeof(double));

    if (!SVD_decomposition->A || !SVD_decomposition->U || !SVD_decomposition->S || !SVD_decomposition->V) {
        fprintf(stderr, "Memory allocation failed for matrices in SVD structure.\n");
        SVD_free(SVD_decomposition);
        return NULL;
    }

    memcpy(SVD_decomposition->A, A, rows * columns * sizeof(double));

    return SVD_decomposition;

}

/*

  #include "LinearAlgebraBasics.h"
#include <math.h> // Pour sqrt()

SVD *SVD_decomposition(double *A, int rows, int columns) {
    // Vérification des arguments
    if (!A || rows <= 0 || columns <= 0) {
        fprintf(stderr, "Error: Invalid input to SVD_decomposition.\n");
        return NULL;
    }

    // Création de la structure SVD
    SVD *SVD_decomposition = create_SVD(A, rows, columns);
    if (!SVD_decomposition) {
        fprintf(stderr, "Error: Failed to allocate memory for SVD structure.\n");
        return NULL;
    }

    // Étape 1 : Calcul de A^T (transpose de A)
    double *AT = malloc(columns * rows * sizeof(double));
    if (!AT) {
        fprintf(stderr, "Error: Memory allocation failed for AT.\n");
        SVD_free(SVD_decomposition);
        return NULL;
    }
    Matrix_transpose(A, rows, columns, AT);

    // Étape 2 : Calcul de A^T * A
    double *ATA = calloc(columns * columns, sizeof(double));
    if (!ATA) {
        fprintf(stderr, "Error: Memory allocation failed for ATA.\n");
        free(AT);
        SVD_free(SVD_decomposition);
        return NULL;
    }

    for (int i = 0; i < columns; i++) {
        for (int j = 0; j < columns; j++) {
            for (int k = 0; k < rows; k++) {
                ATA[i * columns + j] += AT[k * columns + i] * A[k * columns + j];
            }
        }
    }

    // Étape 3 : Calcul des valeurs propres et vecteurs propres de ATA
    double *eigenvalues = calloc(columns, sizeof(double));
    double *eigenvectors = calloc(columns * columns, sizeof(double));
    if (!eigenvalues || !eigenvectors) {
        fprintf(stderr, "Error: Memory allocation failed for eigenvalues or eigenvectors.\n");
        free(AT);
        free(ATA);
        free(eigenvalues);
        free(eigenvectors);
        SVD_free(SVD_decomposition);
        return NULL;
    }

    Matrix_eigenvalues(ATA, columns, eigenvalues, eigenvectors);

    // Remplissage de V et des valeurs singulières (S)
    memcpy(SVD_decomposition->V, eigenvectors, columns * columns * sizeof(double));
    for (int i = 0; i < columns; i++) {
        SVD_decomposition->S[i] = sqrt(eigenvalues[i]);
    }

    // Étape 4 : Calcul de U = A * V * Σ^(-1)
    for (int i = 0; i < rows; i++) {
        for (int j = 0; j < columns; j++) {
            SVD_decomposition->U[i * columns + j] = 0.0;
            for (int k = 0; k < columns; k++) {
                if (SVD_decomposition->S[k] > 1e-10) { // Éviter la division par zéro
                    SVD_decomposition->U[i * columns + j] += A[i * columns + k] *
                                                             SVD_decomposition->V[k * columns + j] /
                                                             SVD_decomposition->S[j];
                }
            }
        }
    }

    // Nettoyage
    free(AT);
    free(ATA);
    free(eigenvalues);
    free(eigenvectors);

    return SVD_decomposition;
}
*/ 

void SVD_free(SVD *SVD_decomposition) {

    if (!SVD_decomposition) return;

    free(SVD_decomposition->A);
    free(SVD_decomposition->U);
    free(SVD_decomposition->S);
    free(SVD_decomposition->V);

    free(SVD_decomposition);
    
}

