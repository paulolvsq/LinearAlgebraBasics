#include "LinearAlgebraBasics.h"

int main() {

    double X[] = {2.0, 3.0, 1.0};
    double Y[] = {1.0, 2.0, 3.0};

    int dimension = 3;

    double *result_vectors_addition = vectors_addition(X, Y, dimension);
    
    printf("##################################### TEST 1 #####################################\n");
    
    for (int i = 0; i < dimension; i++)
	printf("%lf\t", result_vectors_addition[i]);

    printf("\n");

    double *result_scalar_product = scalar_product(X, Y, dimension);

    printf("##################################### TEST 2 #####################################\n");
    
    for (int i = 0; i < dimension; i++)
	printf("%lf\t", result_scalar_product[i]);

    printf("\n");

    double *result_vector_product = vector_product(X, Y);
    
    printf("##################################### TEST 3 #####################################\n");
    
    for (int i = 0; i < dimension; i++)
	printf("%lf\t", result_vector_product[i]);

    printf("\n");

    double norm = vector_norm(X, 3);
    
    printf("##################################### TEST 4 #####################################\n");

    printf("Norm = %lf\n", norm);

    return 0;

}


