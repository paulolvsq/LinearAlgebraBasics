#include "LinearAlgebraBasics.h"

int main() {

    int rows = 6000;
    int columns = 6000;

    double *A = generate_matrix_double(rows, columns);

    printf("##################################### TEST LU DECOMPOSITION SEQUENTIAL #####################################\n");

    time_t begin = time(NULL);
    
    LU *LU_sequential = LU_decomposition(A, rows, columns);
    
    time_t end = time(NULL);

    free(LU_sequential->L);
    free(LU_sequential->U);
    free(LU_sequential);

    printf("Elapsed time : %ld seconds.\n", (end - begin));

    printf("##################################### TEST LU DECOMPOSITION PARALLEL #####################################\n");

    begin = time(NULL);
    
    LU *LU_parallel = LU_decomposition_parallel(A, rows, columns);

    end = time(NULL);

    printf("Elapsed time : %ld seconds.\n", (end - begin));
    
    free(LU_parallel->L);
    free(LU_parallel->U);
    free(LU_parallel);  
    
    return 0;

}
    

    
