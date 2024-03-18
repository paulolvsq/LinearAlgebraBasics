#include "LinearAlgebraBasics.h"

int main() {

    printf("##################################### TEST LU DECOMPOSITION 1 #####################################\n");
    
    int rows = 3, columns = 3;
    double A[] = {1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0};

    LU *LU_test = LU_decomposition(A, rows, columns);

    printf("L =\n");

    for (int i = 0; i < rows; i++) {
	for (int j = 0; j < columns; j++) {
	    printf("%lf\t", LU_test->L[i * columns + j]);
	}
	printf("\n");
    }

    printf("U =\n");

    for (int i = 0; i < rows; i++) {
	for (int j = 0; j < columns; j++) {
	    printf("%lf\t", LU_test->U[i * columns + j]);
	}
	printf("\n");
    }

    free(LU_test->L);
    free(LU_test->U);

    free(LU_test);
    
    printf("##################################### TEST LU DECOMPOSITION 2 #####################################\n");
    
    double B[] = {3.0, 2.0, 4.0, 2.0, 0.0, 2.0, 4.0, 2.0, 3.0};

    LU *LU_test2 = LU_decomposition(B, rows, columns);

    printf("L =\n");

    for (int i = 0; i < rows; i++) {
	for (int j = 0; j < columns; j++) {
	    printf("%lf\t", LU_test2->L[i * columns + j]);
	}
	printf("\n");
    }

    printf("U =\n");

    for (int i = 0; i < rows; i++) {
	for (int j = 0; j < columns; j++) {
	    printf("%lf\t", LU_test2->U[i * columns + j]);
	}
	printf("\n");
    }

    free(LU_test2->L);
    free(LU_test2->U);

    free(LU_test2);
    
    return 0;

}
