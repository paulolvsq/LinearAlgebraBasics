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
    SVD_decomposition->U = calloc(rows * min, sizeof(double));
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

SVD *SVD_decomposition(double *A, int rows, int columns) {
    
    if (!A || rows <= 0 || columns <= 0) {
        fprintf(stderr, "Error: Invalid input to SVD_decomposition.\n");
        return NULL;
    }

    SVD *SVD_decomposition = create_SVD(A, rows, columns);
    
    if (!SVD_decomposition) {
        fprintf(stderr, "Error: Failed to allocate memory for SVD structure.\n");
        return NULL;
    }

    double *AT = matrix_transpose(A, rows, columns);

    if (!AT) {
        fprintf(stderr, "Error: Memory allocation failed for AT.\n");
        SVD_free(SVD_decomposition);
        return NULL;
    }

    printf("AT = \n");

    for (int i = 0; i < columns; i++) {
	for (int j = 0; j < rows; j++) {
	    printf("%lf\t", AT[i * rows + j]);
	}
	printf("\n");
    }
    
    double *ATA = sequential_matrix_product(AT, columns, rows, A, rows, columns); 

    if (!ATA) {
        fprintf(stderr, "Error: Matrix product failed in SVD for ATA.\n");
        free(AT);
        SVD_free(SVD_decomposition);
        return NULL;
    }

    printf("ATA = \n");

    for (int i = 0; i < columns; i++) {
	for (int j = 0; j < columns; j++) {
	    printf("%lf\t", ATA[i * columns + j]);
	}
	printf("\n");
    }
    
    double *eigenvalues = matrix_eigenvalues(ATA, columns, columns, 1000, 1E-20);

    double *eigenvectors = matrix_eigenvectors(ATA, columns, eigenvalues, 1000, 1E-20);

    if (!eigenvalues || !eigenvectors) {
        fprintf(stderr, "Error: Memory allocation failed for eigenvalues or eigenvectors.\n");
        free(AT);
        free(ATA);
        free(eigenvalues);
        free(eigenvectors);
        SVD_free(SVD_decomposition);
        return NULL;
    }

    memcpy(SVD_decomposition->V, eigenvectors, columns * columns * sizeof(double));

    int min = (rows < columns) ? rows : columns;

    for (int i = 0; i < min; i++) {
        SVD_decomposition->S[i] = sqrt(fabs(eigenvalues[i]));
    }

    double epsilon = 1e-12;
    
    for (int j = 0; j < min; j++) {
        if (SVD_decomposition->S[j] > epsilon) {
            for (int i = 0; i < rows; i++) {
                double sum = 0.0;
                for (int k = 0; k < columns; k++) {
                    sum += A[i * columns + k] * SVD_decomposition->V[k * columns + j];
                }
                SVD_decomposition->U[i * min + j] = sum / SVD_decomposition->S[j];
            }
        } else {
            for (int i = 0; i < rows; i++) {
                SVD_decomposition->U[i * min + j] = 0.0;
            }
        }
    }

    free(AT);
    free(ATA);
    free(eigenvalues);
    free(eigenvectors);

    return SVD_decomposition;
}

void SVD_free(SVD *SVD_decomposition) {

    if (!SVD_decomposition) return;

    free(SVD_decomposition->A);
    free(SVD_decomposition->U);
    free(SVD_decomposition->S);
    free(SVD_decomposition->V);

    free(SVD_decomposition);
    
}

