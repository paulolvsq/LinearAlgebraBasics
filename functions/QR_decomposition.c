#include "LinearAlgebraBasics.h"
#include <string.h>

QR *create_QR(double *A, int rows, int columns) {

    QR *QR_decomposition = malloc(sizeof(QR));

    // A (rows * columns) ||| Q (rows * columns) ||| R (columns * columns)
    // A (m * n) ||| Q (m * n) ||| R (n * n)
    // m > n with rank(A) = n
    // to be coherent with memory access in QR_decomposition functions
    // mathematically, Q is (rows * columns) and R is (columns * columns)

    QR_decomposition->rows = rows;
    QR_decomposition->columns = columns;

    QR_decomposition->A = malloc(rows * columns * sizeof(double));

    memcpy(QR_decomposition->A, A, rows * columns * sizeof(double));

    QR_decomposition->Q = malloc(rows * columns * sizeof(double));
    QR_decomposition->R = malloc(columns * columns * sizeof(double));

    return QR_decomposition;

}

QR *QR_decomposition(double *A, int rows, int columns) {
    
    QR *QR_decomposition = create_QR(A, rows, columns);

    // ASSERTION : THE USER WILL ALWAYS GIVE MATRICES OF THE RIGHT SIZE AS INPUT
    // MODIFIED GRAM SCHMIDT METHOD

    double val;
    
    for (int k = 0; k < columns; k++) {
	double s = 0.0;
	for (int j = 0; j < rows; j++) {
	    val = QR_decomposition->A[j * columns + k];
	    s += val * val;
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

    double epsilon = 1e-10; // Pour éviter les divisions par zéro

    // Boucle sur les colonnes pour effectuer la décomposition QR
    for (int k = 0; k < columns; k++) {
        double s = 0.0;

        // Étape 1 : Calcul de la norme ||A[:,k]||^2 en parallèle
#pragma omp parallel for reduction(+:s) schedule(static)
        for (int j = 0; j < rows; j++) {
            double val = QR_decomposition->A[j * columns + k];
            s += val * val;
        }
        QR_decomposition->R[k * columns + k] = sqrt(s);

        // Étape 2 : Normalisation pour obtenir Q[:,k]
#pragma omp parallel for schedule(static)
        for (int j = 0; j < rows; j++) {
            if (QR_decomposition->R[k * columns + k] > epsilon) {
                QR_decomposition->Q[j * columns + k] =
                    QR_decomposition->A[j * columns + k] / QR_decomposition->R[k * columns + k];
            } else {
                QR_decomposition->Q[j * columns + k] = 0.0;
            }
        }

        // Étape 3 : Calcul des coefficients R[k,i]
        for (int i = k + 1; i < columns; i++) {
            s = 0.0;
#pragma omp parallel for reduction(+:s) schedule(static)
            for (int j = 0; j < rows; j++) {
                s += QR_decomposition->A[j * columns + i] *
		    QR_decomposition->Q[j * columns + k];
            }
            QR_decomposition->R[k * columns + i] = s;

            // Étape 4 : Mise à jour de A[:,i]
#pragma omp parallel for schedule(static)
            for (int j = 0; j < rows; j++) {
                QR_decomposition->A[j * columns + i] -=
                    QR_decomposition->R[k * columns + i] *
                    QR_decomposition->Q[j * columns + k];
            }
        }
    }

    return QR_decomposition;
}

/* QR *QR_decomposition_parallel(double *A, int rows, int columns) { */

/*     QR *QR_decomposition = create_QR(A, rows, columns); */
    
/*     // ASSERTION : THE USER WILL ALWAYS GIVE MATRICES OF THE RIGHT SIZE AS INPUT */
/*     // MODIFIED GRAM SCHMIDT METHOD PARALLEL VERSION */

/*     double val; */
    
/*     for (int k = 0; k < columns; k++) { */
/* 	double s = 0.0; */
/* #pragma omp parallel for reduction(+:s) schedule(static) */
/* 	for (int j = 0; j < rows; j++) { */
/* 	    val = QR_decomposition->A[j * columns + k]; */
/* 	    s += val * val; */
/* 	} */
/* 	QR_decomposition->R[k * columns + k] = sqrt(s); */
/* #pragma omp parallel for schedule(static) */
/* 	for (int j = 0; j < rows; j++) { */
/* 	    QR_decomposition->Q[j * columns + k] = QR_decomposition->A[j * columns + k] / QR_decomposition->R[k * columns + k]; */
/* 	} */
/* 	for (int i = k + 1; i < columns; i++) { */
/* 	    s = 0.0; */
/* #pragma omp parallel for reduction(+:s) schedule(static) */
/* 	    for (int j = 0; j < rows; j++) { */
/* 		s += QR_decomposition->A[j * columns + i] * QR_decomposition->Q[j * columns + k]; */
/* 	    } */
/* 	    QR_decomposition->R[k * columns + i] = s; */
/* #pragma omp parallel for schedule(static) */
/* 	    for (int j = 0; j < rows; j++) { */
/* 		QR_decomposition->A[j * columns + i] = QR_decomposition->A[j * columns + i] - QR_decomposition->R[k * columns + i] * QR_decomposition->Q[j * columns + k]; */
/* 	    } */
/* 	} */
/*     } */

/*     return QR_decomposition; */

/* } */

void QR_free(QR *QR_decomposition) {

    free(QR_decomposition->A);
    free(QR_decomposition->Q);
    free(QR_decomposition->R);

    free(QR_decomposition);

}
