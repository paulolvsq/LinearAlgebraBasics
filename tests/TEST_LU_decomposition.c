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

    LU_free(LU_test);
    
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

    LU_free(LU_test2);
    
    printf("##################################### TEST LU DECOMPOSITION PARALLEL 1 #####################################\n");

    LU *LU_test_parallel = LU_decomposition_parallel(A, rows, columns);

    printf("L =\n");

    for (int i = 0; i < rows; i++) {
	for (int j = 0; j < columns; j++) {
	    printf("%lf\t", LU_test_parallel->L[i * columns + j]);
	}
	printf("\n");
    }

    printf("U =\n");

    for (int i = 0; i < rows; i++) {
	for (int j = 0; j < columns; j++) {
	    printf("%lf\t", LU_test_parallel->U[i * columns + j]);
	}
	printf("\n");
    }

    LU_free(LU_test_parallel);
    
    printf("##################################### TEST LU DECOMPOSITION PARALLEL 2 #####################################\n");

    LU *LU_test_parallel2 = LU_decomposition_parallel(B, rows, columns);

    printf("L =\n");

    for (int i = 0; i < rows; i++) {
	for (int j = 0; j < columns; j++) {
	    printf("%lf\t", LU_test_parallel2->L[i * columns + j]);
	}
	printf("\n");
    }

    printf("U =\n");

    for (int i = 0; i < rows; i++) {
	for (int j = 0; j < columns; j++) {
	    printf("%lf\t", LU_test_parallel2->U[i * columns + j]);
	}
	printf("\n");
    }

    LU_free(LU_test_parallel2);
    
    return 0;

}
