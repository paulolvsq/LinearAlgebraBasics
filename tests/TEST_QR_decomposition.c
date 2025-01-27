#include "LinearAlgebraBasics.h"

int main() {

    printf("##################################### TEST QR DECOMPOSITION #####################################\n");
    
    int rows = 3, columns = 3;
    double A[] = {12.0, -51.0, 4.0, 6.0, 167.0, -68.0, -4.0, 24.0, -41.0};

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
    
    printf("##################################### TEST QR DECOMPOSITION PARALLEL #####################################\n");
    
    QR *QR_parallel_test = QR_decomposition_parallel(A, rows, columns);

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

    QR_free(QR_test);
    QR_free(QR_parallel_test);
    
    return 0;

}
