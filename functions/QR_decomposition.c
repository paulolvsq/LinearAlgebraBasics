#include "LinearAlgebraBasics.h"

QR *create_QR(double *A, int rows, int columns) {

    QR *QR_decomposition = malloc(sizeof(QR));

    // A (rows * columns) ||| Q (rows * rows) ||| R (rows * columns)

    QR_decomposition->rows = rows;
    QR_decomposition->columns = columns;
    
    QR_decomposition->A = A;
    QR_decomposition->Q = malloc(rows * rows * sizeof(double));
    QR_decomposition->R = malloc(columns * rows * sizeof(double));

    return QR_decomposition;

}
