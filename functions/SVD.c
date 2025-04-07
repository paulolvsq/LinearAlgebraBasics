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

void SVD_free(SVD *SVD_decomposition) {

    if (!SVD_decomposition) return;

    free(SVD_decomposition->A);
    free(SVD_decomposition->U);
    free(SVD_decomposition->S);
    free(SVD_decomposition->V);

    free(SVD_decomposition);
    
}

