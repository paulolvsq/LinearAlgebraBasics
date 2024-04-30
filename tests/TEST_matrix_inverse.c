#include "LinearAlgebraBasics.h"

int main() {

    double A[] = {1.0, 1.0, 1.0, 3.0, 1.0, -3.0, 1.0, -2.0, -5.0};
    double b[] = {1.0, 5.0, 10.0};
    
    int d = 3;

    double *res = solve_LU_system(A, b, d);
        
    printf("##################################### TEST LU SOLVE #####################################\n");

    for (int i = 0; i < d; i++)
	printf("%lf\t", res[i]);
    
    printf("\n");

    printf("##################################### TEST MATRIX INVERSE #####################################\n");

    double *inverse = matrix_inverse(A, d);
    
    for (int i = 0; i < d; i++) {
	for (int j = 0; j < d; j++) {
	    printf("%lf\t", inverse[i * d + j]);
	}
	printf("\n");
    }

    return 0;

}
