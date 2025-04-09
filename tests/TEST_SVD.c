#include "LinearAlgebraBasics.h"

int main() {

    double A[] = {3.0, 2.0, 2.0,
	2.0, 3.0, -2.0};

    int rows = 2;
    int columns = 3;

    SVD *SVD_test = SVD_decomposition(A, rows, columns);

    printf("################################### TEST SVD ###################################\n");

    int min = (rows < columns) ? rows : columns;

    printf("U = \n");

    for (int i = 0; i < rows; i++) {
	for (int j = 0; j < min; j++) {
	    printf("%lf\t", SVD_test->U[i * min + j]);
	}
	printf("\n");
    }

    printf("S = \n");

    for (int i = 0; i < min; i++) {
	printf("%lf\t", SVD_test->S[i]);
    }

    printf("V = \n");

    for (int i = 0; i < columns; i++) {
	for (int j = 0; j < columns; j++) {
	    printf("%lf\t", SVD_test->V[i * columns + j]);
	}
	printf("\n");
    }

    SVD_free(SVD_test);
    
    return 0;

}
