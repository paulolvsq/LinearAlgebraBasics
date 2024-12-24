#include "LinearAlgebraBasics.h"

int main() {

    int P_rows = 3;
    int P_columns = 3;

    int Q_rows = 3;
    int Q_columns = 2;

    double P[] = {1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0};
    double Q[] = {1.0, 2.0, 3.0, 4.0, 5.0, 6.0};
    
    double *matrix = parallel_matrix_product(P, P_rows, P_columns, Q, Q_rows, Q_columns);

    printf("##################################### TEST 1 #####################################\n");
    
    for (int i = 0; i < P_rows; i++) {
	for (int j = 0; j < Q_columns; j++) {
	    printf("%lf\t", matrix[i * Q_columns + j]);
	}
	printf("\n");
    }

    free(matrix);
    
    P_rows = 4;
    P_columns = 5;

    Q_rows = 5;
    Q_columns = 2;

    double P2[] = {1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0, 10.0, 11.0, 12.0, 13.0, 14.0, 15.0, 16.0, 17.0, 18.0, 19.0, 20.0};
    double Q2[] = {1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0, 10.0};
    
    double *matrix2 = parallel_matrix_product(P2, P_rows, P_columns, Q2, Q_rows, Q_columns);

    printf("##################################### TEST 2 #####################################\n");
    
    for (int i = 0; i < P_rows; i++) {
	for (int j = 0; j < Q_columns; j++) {
	    printf("%lf\t", matrix2[i * Q_columns + j]);
	}
	printf("\n");
    }

    free(matrix2);
    
    P_rows = 4;
    P_columns = 5;

    int dimension = 5;
    
    double P3[] = {1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0, 10.0, 11.0, 12.0, 13.0, 14.0, 15.0, 16.0, 17.0, 18.0, 19.0, 20.0};
    double *Q3 = generate_identity_matrix(dimension);

    printf("##################################### TEST 3 #####################################\n");

    // Matrix of size (4,5) x (5,5) id => matrix of size (4,5)
    
    double *matrix3 = parallel_matrix_product(P3, P_rows, P_columns, Q3, dimension, dimension);
    
    for (int i = 0; i < P_rows; i++) {
	for (int j = 0; j < P_columns; j++) {
 	    printf("%lf\t", matrix3[i * dimension + j]);
	}
	printf("\n");
    }

    free(Q3);
    free(matrix3);
    
    return 0;

}
