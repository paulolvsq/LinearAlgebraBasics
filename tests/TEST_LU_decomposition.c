#include "LinearAlgebraBasics.h"

int main() {

    printf("##################################### TEST LU DECOMPOSITION #####################################\n");
    
    int dimension = 3;
    double A[] = {1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0};

    LU *LU_test = LU_decomposition(A, dimension);

    printf("L =\n");

    for (int i = 0; i < dimension; i++) {
	for (int j = 0; j < dimension; j++) {
	    printf("%lf\t", LU_test->L[i * dimension + j]);
	}
	printf("\n");
    }

    printf("U =\n");

    for (int i = 0; i < dimension; i++) {
	for (int j = 0; j < dimension; j++) {
	    printf("%lf\t", LU_test->U[i * dimension + j]);
	}
	printf("\n");
    }

    free(LU_test->L);
    free(LU_test->U);
    
    return 0;

}
