#include "LinearAlgebraBasics.h"
#define DOUBLE 100.0

int main() {

    int A_rows = 40000;
    int A_columns = 40000;
    int dimension = 40000;

    /* I can't test with more, otherwise I exceed the memory size of an integer when I allocate the 
       matrix with a (rows x columns) size */
    
    double *A = generate_matrix_double(A_rows, A_columns);
    double *X = malloc(dimension * sizeof(double));

    srand(time(NULL));

    for (int i = 0; i < dimension; i++)
	X[i] = (double) (rand()) / RAND_MAX * DOUBLE;

    printf("##################################### TEST SEQUENTIAL VECTOR MATRIX PRODUCT #####################################\n");
    
    time_t begin = time(NULL);
    
    double *vector_sequential = sequential_vector_matrix_product(A, A_rows, A_columns, X, dimension);

    time_t end = time(NULL);

    printf("Elapsed time : %ld seconds.\n", (end - begin));
    
    free(vector_sequential);
    
    printf("##################################### TEST PARALLEL VECTOR MATRIX PRODUCT #####################################\n");

    begin = time(NULL);
    
    double *vector_parallel = parallel_vector_matrix_product(A, A_rows, A_columns, X, dimension);

    end = time(NULL);

    printf("Elapsed time : %ld seconds.\n", (end - begin));
    
    free(vector_parallel);

    free(X);
    free(A);
    
    return 0;

}
