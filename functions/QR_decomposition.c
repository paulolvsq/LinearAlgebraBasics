#include "LinearAlgebraBasics.h"
#include <string.h>

QR *create_QR(double *A, int rows, int columns) {

    QR *QR_decomposition = malloc(sizeof(QR));

    // A (rows * columns) ||| Q (rows * rows) ||| R (rows * columns)
    // A (m * n) ||| Q (m * m) ||| R (m * n)
    // m > n with rank(A) = n

    QR_decomposition->rows = rows;
    QR_decomposition->columns = columns;

    QR_decomposition->A = malloc(rows * columns * sizeof(double));

    memcpy(QR_decomposition->A, A, rows * columns * sizeof(double));

    QR_decomposition->Q = malloc(rows * rows * sizeof(double));
    QR_decomposition->R = malloc(columns * rows * sizeof(double));

    return QR_decomposition;

}

QR *QR_decomposition(double *A, int rows, int columns) {

    QR *QR_decomposition = create_QR(A, rows, columns);

    // ASSERTION : THE USER WILL ALWAYS GIVE MATRICES OF THE RIGHT SIZE AS INPUT
    // MODIFIED GRAM SCHMIDT METHOD

    for (int k = 0; k < columns; k++) {
	double s = 0.0;
	for (int j = 0; j < rows; j++) {
	    s += pow(QR_decomposition->A[j * columns + k], 2);
	}
	QR_decomposition->R[k * columns + k] = sqrt(s);
	for (int j = 0; j < rows; j++) {
	    QR_decomposition->Q[j * columns + k] = QR_decomposition->A[j * columns + k] / QR_decomposition->R[k * columns + k];
	}
	for (int i = k + 1; i < columns; i++) {
	    s = 0.0;
	    for (int j = 0; j < rows; j++) {
		s += QR_decomposition->A[j * columns + i] * QR_decomposition->Q[j * columns + k];
	    }
	    QR_decomposition->R[k * columns + i] = s;
	    for (int j = 0; j < rows; j++) {
		QR_decomposition->A[j * columns + i] = QR_decomposition->A[j * columns + i] - QR_decomposition->R[k * columns + i] * QR_decomposition->Q[j * columns + k];
	    }
	}
    }

    return QR_decomposition;
    
}

QR *QR_decomposition_parallel(double *A, int rows, int columns) {

    QR *QR_decomposition = create_QR(A, rows, columns);
    
    // ASSERTION : THE USER WILL ALWAYS GIVE MATRICES OF THE RIGHT SIZE AS INPUT
    // MODIFIED GRAM SCHMIDT METHOD PARALLEL VERSION

    for (int k = 0; k < columns; k++) {
	double s = 0.0;
#pragma omp parallel for reduction(+:s)
	for (int j = 0; j < rows; j++) {
	    s += pow(QR_decomposition->A[j * columns + k], 2);
	}
	QR_decomposition->R[k * columns + k] = sqrt(s);
#pragma omp parallel for
	for (int j = 0; j < rows; j++) {
	    QR_decomposition->Q[j * columns + k] = QR_decomposition->A[j * columns + k] / QR_decomposition->R[k * columns + k];
	}
	for (int i = k + 1; i < columns; i++) {
	    s = 0.0;
#pragma omp parallel for reduction(+:s)
	    for (int j = 0; j < rows; j++) {
		s += QR_decomposition->A[j * columns + i] * QR_decomposition->Q[j * columns + k];
	    }
	    QR_decomposition->R[k * columns + i] = s;
#pragma omp parallel for
	    for (int j = 0; j < rows; j++) {
		QR_decomposition->A[j * columns + i] = QR_decomposition->A[j * columns + i] - QR_decomposition->R[k * columns + i] * QR_decomposition->Q[j * columns + k];
	    }
	}
    }

    return QR_decomposition;

}

void QR_free(QR *QR_decomposition) {

    free(QR_decomposition->A);
    free(QR_decomposition->Q);
    free(QR_decomposition->R);

    free(QR_decomposition);

}
