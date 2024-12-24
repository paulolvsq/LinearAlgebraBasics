#include "LinearAlgebraBasics.h"
#include <string.h>

int main() {

    int rows = 3000;
    int columns = 3000;

    double *A = generate_matrix_double(rows, columns);
    double *B = malloc((rows * columns) * sizeof(double));
    memcpy(B, A, rows * columns * sizeof(double));

    printf("##################################### TEST QR DECOMPOSITION SEQUENTIAL #####################################\n");

    time_t begin = time(NULL);
    
    QR *QR_sequential = QR_decomposition(A, rows, columns);
    
    time_t end = time(NULL);

    printf("Elapsed time : %ld seconds.\n", (end - begin));
    
    printf("##################################### TEST QR DECOMPOSITION PARALLEL #####################################\n");

    begin = time(NULL);
    
    QR *QR_parallel = QR_decomposition_parallel(B, rows, columns);

    end = time(NULL);

    printf("Elapsed time : %ld seconds.\n", (end - begin));
    
    free(B);
    
    free(QR_parallel->Q);
    free(QR_parallel->R);
    free(QR_parallel);  

    free(A);
    
    free(QR_sequential->Q);
    free(QR_sequential->R);
    free(QR_sequential);

    
    return 0;

}
