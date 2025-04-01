#include "LinearAlgebraBasics.h"

int main() {

    double A[] = {1.0, 1.0, 1.0, 1.0,
	1.0, 5.0, 5.0, 5.0,
	1.0, 5.0, 14.0, 14.0,
	1.0, 5.0, 14.0, 15.0};

    int size = 4;

    Cholesky *Cholesky_test = Cholesky_decomposition(A, size);

    printf("############################# TEST CHOLESKY #############################\n");

    printf("L = \n");
    
    for (int i = 0; i < size; i++) {
	for (int j = 0; j < size; j++) {
	    printf("%lf\t", Cholesky_test->L[i * size + j]);
	}
	printf("\n");
    }

    printf("L_t = \n");
    
    for (int i = 0; i < size; i++) {
	for (int j = 0; j < size; j++) {
	    printf("%lf\t", Cholesky_test->L_t[i * size + j]);
	}
	printf("\n");
    }

    free_Cholesky(Cholesky_test);

    return 0;
    
}

    
