#include "LinearAlgebraBasics.h"

int main() {

    int P_rows = 3;
    int P_columns = 3;

    double P[] = {1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0};

    int scalar = 3;

    double *result_matrices_addition = matrices_addition(P, P, P_rows, P_columns);
    
    printf("##################################### TEST 1 #####################################\n");
    
    for (int i = 0; i < P_rows; i++) {
	for (int j = 0; j < P_columns; j++) {
	    printf("%lf\t", result_matrices_addition[i * P_columns + j]);
	}
	printf("\n");
    }
    printf("\n");

    double *result_scalar_multiplication = matrix_scalar_multiplication(P, P_rows, P_columns, scalar);

    printf("##################################### TEST 2 #####################################\n");
    
    for (int i = 0; i < P_rows; i++) {
	for (int j = 0; j < P_columns; j++) {
	    printf("%lf\t", result_scalar_multiplication[i * P_columns + j]);
	}
	printf("\n");
    }
    printf("\n");

    double *result_matrix_transpose = matrix_transpose(P, P_rows, P_columns);
    
    printf("##################################### TEST 3 #####################################\n");
        
    for (int i = 0; i < P_rows; i++) {
	for (int j = 0; j < P_columns; j++) {
	    printf("%lf\t", result_matrix_transpose[i * P_columns + j]);
	}
	printf("\n");
    }
    printf("\n");

    int dimension = 10;
    double *identity = generate_identity_matrix(dimension);
    double trace = matrix_trace(identity, dimension);

    printf("##################################### TEST 4 #####################################\n");
    
    printf("Trace = %lf\n", trace);

    return 0;

}
