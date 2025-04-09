#include "LinearAlgebraBasics.h"
#include <string.h>

/**
 * @brief Allocates and initializes a QR decomposition structure.
 *
 * This function creates a QR decomposition structure for a given matrix A
 * and initializes its components (Q, R, and A).
 *
 * @param A Pointer to the input matrix (size: rows x columns).
 * @param rows Number of rows in the matrix (must be positive).
 * @param columns Number of columns in the matrix (must be positive).
 *
 * @return Pointer to the QR structure containing Q, R, and A matrices on success,
 *         or NULL on failure due to invalid dimensions or memory allocation errors.
 */

QR *create_QR(double *A, int rows, int columns) {

    if (rows <= 0 || columns <= 0) {
        fprintf(stderr, "Error: Invalid dimensions for QR decomposition (rows=%d, columns=%d). Both must be strictly positive.\n", rows, columns);
        return NULL;
    }

    if (!A) {
        fprintf(stderr, "Error: Null pointer detected for input matrix in create_QR.\n");
        return NULL;
    }

    QR *QR_decomposition = malloc(sizeof(QR));

    if (!QR_decomposition) {
        fprintf(stderr, "Error: Memory allocation failed for QR structure.\n");
        return NULL;
    }

    // A (rows * columns) ||| Q (rows * columns) ||| R (columns * columns)
    // A (m * n) ||| Q (m * n) ||| R (n * n)
    // m > n with rank(A) = n
    // to be coherent with memory access in QR_decomposition functions
    // mathematically, Q is (rows * columns) and R is (columns * columns)

    QR_decomposition->rows = rows;
    QR_decomposition->columns = columns;

    QR_decomposition->A = malloc(rows * columns * sizeof(double));
    QR_decomposition->Q = malloc(rows * columns * sizeof(double));
    QR_decomposition->R = calloc(columns * columns, sizeof(double));

    if (!QR_decomposition->A || !QR_decomposition->Q || !QR_decomposition->R) {
        fprintf(stderr, "Error: Memory allocation failed for matrices in create_QR.\n");
        free(QR_decomposition->A);
        free(QR_decomposition->Q);
        free(QR_decomposition->R);
        free(QR_decomposition);
        return NULL;
    }

    memcpy(QR_decomposition->A, A, rows * columns * sizeof(double));
    
    return QR_decomposition;

}

/**
 * @brief Performs QR decomposition on a given matrix using the Modified Gram-Schmidt method.
 *
 * This function decomposes a matrix A into an orthogonal matrix Q
 * and an upper triangular matrix R such that A = Q * R.
 *
 * @param A Pointer to the input matrix (size: rows x columns).
 * @param rows Number of rows in the matrix (must be positive).
 * @param columns Number of columns in the matrix (must be positive).
 *
 * @return Pointer to the QR structure containing Q and R matrices on success,
 *         or NULL on failure due to invalid dimensions, singular columns, or memory allocation errors.
 */

QR *QR_decomposition(double *A, int rows, int columns) {
    
    QR *QR_decomposition = create_QR(A, rows, columns);

    if (!QR_decomposition) {
        fprintf(stderr, "Error: Failed to create QR decomposition structure.\n");
        return NULL;
    }

    // MODIFIED GRAM SCHMIDT METHOD

    double val;
    double epsilon = 1e-10;
    
    for (int k = 0; k < columns; k++) {

	double s = 0.0;

	for (int j = 0; j < rows; j++) {
	    val = QR_decomposition->A[j * columns + k];
	    printf("val = %lf\n", val);
	    s += val * val;
	}
	
	if (s < epsilon) { 
            fprintf(stderr, "Error: Column %d is singular or zero during QR decomposition and s = %lf.\n", k, s);
            QR_free(QR_decomposition);
            return NULL;
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

/**
 * @brief Performs parallelized QR decomposition on a given matrix using OpenMP.
 *
 * This function decomposes a matrix A into an orthogonal matrix Q
 * and an upper triangular matrix R such that A = Q * R. The computation is
 * parallelized using OpenMP for improved performance on large matrices.
 *
 * @param A Pointer to the input matrix (size: rows x columns).
 * @param rows Number of rows in the matrix (must be positive).
 * @param columns Number of columns in the matrix (must be positive).
 *
 * @return Pointer to the QR structure containing Q and R matrices on success,
 *         or NULL on failure due to invalid dimensions, singular columns, or memory allocation errors.
 */

QR *QR_decomposition_parallel(double *A, int rows, int columns) {
    
    QR *QR_decomposition = create_QR(A, rows, columns);

    if (!QR_decomposition) {
        fprintf(stderr, "Error: Failed to create parallelized QR decomposition structure.\n");
        return NULL;
    }

    // MODIFIED GRAM SCHMIDT METHOD
    
    double epsilon = 1e-10; 

    for (int k = 0; k < columns; k++) {

	double s = 0.0;

#pragma omp parallel for reduction(+:s) schedule(static)
        for (int j = 0; j < rows; j++) {
	    double val = QR_decomposition->A[j * columns + k];
            s += val * val;
        }
	
	if (s < epsilon) {
            fprintf(stderr, "Error: Singular column detected at column %d during QR decomposition.\n", k);
            QR_free(QR_decomposition);
            return NULL;
        }
	
        QR_decomposition->R[k * columns + k] = sqrt(s);

#pragma omp parallel for schedule(static)
        for (int j = 0; j < rows; j++) {
            if (QR_decomposition->R[k * columns + k] > epsilon) {
                QR_decomposition->Q[j * columns + k] =
                    QR_decomposition->A[j * columns + k] / QR_decomposition->R[k * columns + k];
            } else {
                QR_decomposition->Q[j * columns + k] = 0.0;
            }
        }
	
        for (int i = k + 1; i < columns; i++) {
            s = 0.0;
#pragma omp parallel for reduction(+:s) schedule(static)
            for (int j = 0; j < rows; j++) {
                s += QR_decomposition->A[j * columns + i] * QR_decomposition->Q[j * columns + k];
            }
            QR_decomposition->R[k * columns + i] = s;
#pragma omp parallel for schedule(static)
            for (int j = 0; j < rows; j++) {
                QR_decomposition->A[j * columns + i] -= QR_decomposition->R[k * columns + i] * QR_decomposition->Q[j * columns + k];
            }
        }
	
    }

    return QR_decomposition;
}

/**
 * @brief Frees all memory associated with a QR decomposition structure.
 *
 * This function releases the memory allocated for the matrices Q, R, and A,
 * as well as the QR structure itself.
 *
 * @param QR_decomposition Pointer to the QR structure to free.
 */

void QR_free(QR *QR_decomposition) {

    if (!QR_decomposition) return;

    free(QR_decomposition->A);
    free(QR_decomposition->Q);
    free(QR_decomposition->R);

    free(QR_decomposition);

}
