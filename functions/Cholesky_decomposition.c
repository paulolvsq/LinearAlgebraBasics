#include "LinearAlgebraBasics.h"
#include <string.h>

Cholesky *create_Cholesky(double *A, int size) {
    
    Cholesky *Cholesky_decomposition = malloc(sizeof(Cholesky));

    Cholesky_decomposition->size = size;

    Cholesky_decomposition->A = malloc(size * size * sizeof(double));
    Cholesky_decomposition->L = calloc(size * size, sizeof(double));
    Cholesky_decomposition->L_t = malloc(size * size * sizeof(double));

    memcpy(Cholesky_decomposition->A, A, size * size * sizeof(double));

    return Cholesky_decomposition;

}

Cholesky *Cholesky_decomposition(double *A, int size) {

    Cholesky *Cholesky_decomposition = create_Cholesky(A, size);

    for (int i = 0; i < size; i++) {
	for (int j = 0; j < size; j++) {

	    int sum = 0;

	    if (j == i) {
		for (int k = 0; k < j; k++)
		    sum += Cholesky_decomposition->L[k * 
    
}



    
