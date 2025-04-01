#include "LinearAlgebraBasics.h"

int main() {

    double A[] = {7.0, 3.0, -1.0, 2.0,
	3.0, 8.0, 1.0, -4.0,
	-1.0, 1.0, 4.0, -1.0,
	2.0, -4.0, -1.0, 6.0};
    
    int size = 4;

    LDLT *LDLT_test = LDLT_decomposition(A, size);

    printf("############################# TEST LDLT #############################\n");

    printf("L = \n");
    
    for (int i = 0; i < size; i++) {
	for (int j = 0; j < size; j++) {
	    printf("%lf\t", LDLT_test->L[i * size + j]);
	}
	printf("\n");
    }

    printf("D = \n");

    for (int i = 0; i < size; i++) {
	for (int j = 0; j < size; j++) {
	    printf("%lf\t", LDLT_test->D[i * size + j]);
	}
	printf("\n");
    }
    
    printf("L_t = \n");
    
    for (int i = 0; i < size; i++) {
	for (int j = 0; j < size; j++) {
	    printf("%lf\t", LDLT_test->L_t[i * size + j]);
	}
	printf("\n");
    }

    free_LDLT(LDLT_test);

    return 0;
    
}
