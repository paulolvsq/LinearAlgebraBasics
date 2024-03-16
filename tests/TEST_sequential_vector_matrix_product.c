#include "LinearAlgebraBasics.h"

int main() {

    int A_rows = 3;
    int A_columns = 3;
    int dimension = 3;
    
    double A[] = {1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0};
    double X[] = {1.0, 2.0, 3.0};

    double *vector = sequential_vector_matrix_product(A, A_rows, A_columns, X, dimension);

    printf("##################################### TEST 1 #####################################\n");
    
    for (int i = 0; i < A_rows; i++)
	printf("%lf\t", vector[i]);

    printf("\n");

    free(vector);
    
    A_rows = 5;
    A_columns = 4;
    dimension = 4;
    
    double A2[] = {1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0, 10.0, 11.0, 12.0, 13.0, 14.0, 15.0, 16.0, 17.0, 18.0, 19.0, 20.0};
    double X2[] = {1.0, 2.0, 3.0, 4.0};

    double *vector2 = sequential_vector_matrix_product(A2, A_rows, A_columns, X2, dimension);

    printf("##################################### TEST 2 #####################################\n");
    
    for (int i = 0; i < A_rows; i++)
	printf("%lf\t", vector2[i]);

    free(vector2);
    
    printf("\n");

    dimension = 4;
    
    double *A3 = generate_identity_matrix(dimension);
    double X3[] = {1.0, 2.0, 3.0, 4.0};

    double *vector3 = sequential_vector_matrix_product(A3, dimension, dimension, X3, dimension);

    printf("##################################### TEST 3 #####################################\n");
    
    for (int i = 0; i < dimension; i++)
	printf("%lf\t", vector3[i]);

    printf("\n");

    free(A3);
    free(vector3);

    return 0;

}
