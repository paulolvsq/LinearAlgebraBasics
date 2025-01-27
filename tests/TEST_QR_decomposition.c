#include "LinearAlgebraBasics.h"
#include <string.h>

int main() {

    printf("##################################### TEST QR DECOMPOSITION #####################################\n");
    
    int rows = 3, columns = 3;
    double A[] = {12.0, -51.0, 4.0, 6.0, 167.0, -68.0, -4.0, 24.0, -41.0};

    double *B = malloc((rows * columns) * sizeof(double));
    memcpy(B, A, rows * columns * sizeof(double));
    
    QR *QR_test = QR_decomposition(A, rows, columns);

    printf("Q =\n");

    for (int i = 0; i < rows; i++) {
	for (int j = 0; j < columns; j++) {
	    printf("%lf\t", QR_test->Q[i * columns + j]);
	}
	printf("\n");
    }

    printf("R =\n");

    for (int i = 0; i < rows; i++) {
	for (int j = 0; j < columns; j++) {
	    printf("%lf\t", QR_test->R[i * columns + j]);
	}
	printf("\n");
    }

    free(QR_test->Q);
    free(QR_test->R);
    free(QR_test);
    
    printf("##################################### TEST QR DECOMPOSITION PARALLEL #####################################\n");
    
    QR *QR_parallel_test = QR_decomposition_parallel(B, rows, columns);

    printf("Q =\n");

    for (int i = 0; i < rows; i++) {
	for (int j = 0; j < columns; j++) {
	    printf("%lf\t", QR_parallel_test->Q[i * columns + j]);
	}
	printf("\n");
    }

    printf("R =\n");

    for (int i = 0; i < rows; i++) {
	for (int j = 0; j < columns; j++) {
	    printf("%lf\t", QR_parallel_test->R[i * columns + j]);
	}
	printf("\n");
    }
    
    free(B);
    
    free(QR_parallel_test->Q);
    free(QR_parallel_test->R);
    free(QR_parallel_test);    
    
    return 0;

}
