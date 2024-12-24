#include "LinearAlgebraBasics.h"

int main() {

    int P_rows = 2000;
    int P_columns = 2000;
    int Q_rows = 2000;
    int Q_columns = 2000;

    double *P = generate_matrix_double(P_rows, P_columns);
    double *Q = generate_matrix_double(Q_rows, Q_columns);

    printf("##################################### TEST SEQUENTIAL MATRIX PRODUCT #####################################\n");

    time_t begin = time(NULL);

    double *X_sequential = sequential_matrix_product(P, P_rows, P_columns, Q, Q_rows, Q_columns);

    time_t end = time(NULL);

    printf("Elapsed time : %ld seconds.\n", (end - begin));

    printf("##################################### TEST PARALLEL MATRIX PRODUCT #####################################\n");

    begin = time(NULL);

    double *X_parallel = parallel_matrix_product(P, P_rows, P_columns, Q, Q_rows, Q_columns);

    end = time(NULL);

    printf("Elapsed time : %ld seconds.\n", (end - begin));

    free(P);
    free(Q);
    free(X_sequential);
    free(X_parallel);
    
    return 0;
    
}
