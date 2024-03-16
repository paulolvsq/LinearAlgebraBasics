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

    free(result_matrices_addition);
    
    double *result_scalar_multiplication = matrix_scalar_multiplication(P, P_rows, P_columns, scalar);

    printf("##################################### TEST 2 #####################################\n");
    
    for (int i = 0; i < P_rows; i++) {
	for (int j = 0; j < P_columns; j++) {
	    printf("%lf\t", result_scalar_multiplication[i * P_columns + j]);
	}
	printf("\n");
    }
    printf("\n");

    free(result_scalar_multiplication);

    double *result_matrix_transpose = matrix_transpose(P, P_rows, P_columns);
    
    printf("##################################### TEST 3 #####################################\n");
        
    for (int i = 0; i < P_rows; i++) {
	for (int j = 0; j < P_columns; j++) {
	    printf("%lf\t", result_matrix_transpose[i * P_columns + j]);
	}
	printf("\n");
    }
    printf("\n");

    free(result_matrix_transpose);
    
    int dimension = 10;
    double *identity = generate_identity_matrix(dimension);
    double trace = matrix_trace(identity, dimension);

    printf("##################################### TEST 4 #####################################\n");
    
    printf("Trace = %lf\n", trace);

    free(identity);
    
    printf("##################################### TEST 5 #####################################\n");

    int rows = 3;
    int columns = 2;
    double M[] = {1.0, 2.0, 3.0, 4.0, 5.0, 6.0};

    double norm = matrix_norm(M, rows, columns);

    printf("Norm = %lf\n", norm);
	
    printf("##################################### TEST 6 #####################################\n");
	
    double frobenius = frobenius_norm(M, rows, columns);

    printf("Frobenius = %lf\n", frobenius);

    printf("##################################### TEST 7 #####################################\n");

    double determinant = matrix_determinant(P, P_rows);

    printf("Determinant = %lf\n", determinant);

    printf("##################################### TEST 8 #####################################\n");

    double B[] = {3.0, 2.0, 4.0, 2.0, 0.0, 2.0, 4.0, 2.0, 3.0};

    double determinant2 = matrix_determinant(B, P_rows);

    printf("Determinant = %lf\n", determinant2);
    
    return 0;

}
