#ifndef __LinearAlgebraBasics_
#define __LinearAlgebraBasics_

#include <stdio.h>
#include <stdlib.h>
#include <time.h>

#include <math.h>

/**
 * @brief Represents the LU decomposition of a matrix.
 *
 * This structure stores the components of an LU decomposition, where:
 * - A is the original matrix,
 * - L is the lower triangular matrix,
 * - U is the upper triangular matrix.
 *
 * @struct LU
 * @var LU::rows
 * Number of rows in the matrix.
 * @var LU::columns
 * Number of columns in the matrix.
 * @var LU::A
 * Pointer to the original matrix (size: rows x columns).
 * @var LU::L
 * Pointer to the lower triangular matrix (size: rows x columns).
 * @var LU::U
 * Pointer to the upper triangular matrix (size: rows x columns).
 */

typedef struct LU {

    int rows, columns;

    double *A;
    double *L;
    double *U;

} LU;

/**
 * @brief Represents the QR decomposition of a matrix.
 *
 * This structure stores the components of a QR decomposition, where:
 * - A is the original matrix,
 * - Q is the orthogonal matrix,
 * - R is the upper triangular matrix.
 *
 * @struct QR
 * @var QR::rows
 * Number of rows in the matrix.
 * @var QR::columns
 * Number of columns in the matrix.
 * @var QR::A
 * Pointer to the original matrix (size: rows x columns).
 * @var QR::Q
 * Pointer to the orthogonal matrix (size: rows x columns).
 * @var QR::R
 * Pointer to the upper triangular matrix (size: columns x columns).
 */

typedef struct QR {

    int rows, columns;

    double *A;
    double *Q;
    double *R;

} QR;

/**
 * @brief Represents the Cholesky decomposition of a matrix.
 *
 * This structure stores the components of a Cholesky decomposition, where:
 * - A is the original square matrix,
 * - L is the lower triangular matrix resulting from the decomposition,
 * - Lᵀ is the transpose of the lower triangular matrix.
 *
 * @struct Cholesky
 * @var Cholesky::size
 * Dimension of the square matrix A (size x size).
 * @var Cholesky::A
 * Pointer to the original square matrix (size: size x size).
 * @var Cholesky::L
 * Pointer to the lower triangular matrix resulting from the decomposition (size: size x size).
 * @var Cholesky::L_t
 * Pointer to the transpose of the lower triangular matrix (size: size x size).
 */

typedef struct Cholesky {

    int size;

    double *A;
    double *L;
    double *L_t;

} Cholesky;

/**
 * @brief Represents the LDLT decomposition of a symmetric square matrix.
 *
 * This structure stores the components of an LDLT decomposition, where:
 * - A is the original square matrix,
 * - L is the lower triangular matrix with unit diagonal entries,
 * - D is the diagonal matrix containing the diagonal elements of the decomposition,
 * - Lᵀ is the transpose of the lower triangular matrix.
 *
 * The decomposition satisfies \( A = L \cdot D \cdot L^T \).
 *
 * @struct LDLT
 * @var LDLT::size
 * Dimension of the square matrix A (size x size).
 * @var LDLT::A
 * Pointer to the original square matrix (size: size x size).
 * @var LDLT::L
 * Pointer to the lower triangular matrix resulting from the decomposition (size: size x size).
 * @var LDLT::D
 * Pointer to the diagonal matrix containing diagonal elements (size: size x size).
 * @var LDLT::L_t
 * Pointer to the transpose of the lower triangular matrix (size: size x size).
 */

typedef struct LDLT {

    int size;

    double *A;
    double *L;
    double *D;
    double *L_t;

} LDLT;

typedef struct SVD {

    int rows;
    int columns;

    double *A;
    double *U;
    double *S;
    double *V;

} SVD;

/* SVD.c */

SVD *create_SVD(double *A, int rows, int columns);
SVD *SVD_decomposition(double *A, int rows, int columns);
void SVD_free(SVD *SVD_decomposition);

/* LDLT_decomposition.c */

/**
 * @brief Allocates and initializes an LDLT decomposition structure.
 *
 * This function creates an LDLT decomposition structure for a given matrix A
 * and initializes its components (A, L, D, and Lᵀ). The input matrix A is copied
 * into the structure, and the matrices L, D, and Lᵀ are initialized to zero.
 *
 * @param A Pointer to the input square matrix (size: size x size).
 * @param size Dimension of the square matrix A (must be positive).
 *
 * @return Pointer to the LDLT structure on success, or NULL on failure due to
 *         invalid dimensions or memory allocation errors.
 */

LDLT *create_LDLT(double *A, int size);

/**
 * @brief Performs LDLT decomposition on a symmetric square matrix.
 *
 * This function decomposes a symmetric square matrix A into:
 * - L: a lower triangular matrix with unit diagonal,
 * - D: a diagonal matrix,
 * - Lᵀ: the transpose of the lower triangular matrix.
 *
 * The decomposition satisfies \( A = L \cdot D \cdot L^T \).
 *
 * The function checks if the input matrix is symmetric before performing the decomposition.
 *
 * @param A Pointer to the input square matrix (size: size x size).
 * @param size Dimension of the square matrix A (must be positive).
 *
 * @return Pointer to the LDLT structure containing matrices A, L, D, and Lᵀ on success,
 *         or NULL on failure due to invalid dimensions, non-symmetric matrix,
 *         singularity detection, or memory allocation errors.
 */

LDLT *LDLT_decomposition(double *A, int size);

/**
 * @brief Frees all memory associated with an LDLT decomposition structure.
 *
 * This function releases the memory allocated for the matrices A, L, D, and Lᵀ,
 * as well as the LDLT structure itself.
 *
 * @param LDLT_decomposition Pointer to the LDLT structure to free.
 */

void free_LDLT(LDLT *LDLT_decomposition);

/* Cholesky_decomposition.c */

/**
 * @brief Allocates and initializes a Cholesky decomposition structure.
 *
 * This function creates a Cholesky decomposition structure for a given matrix A
 * and initializes its components (A, L, and Lᵀ). The input matrix A is copied into
 * the structure, and the matrices L and Lᵀ are initialized to zero.
 *
 * @param A Pointer to the input square matrix (size: size x size).
 * @param size Dimension of the square matrix A (must be positive).
 *
 * @return Pointer to the Cholesky structure on success, or NULL on failure due to
 *         invalid dimensions or memory allocation errors.
 */

Cholesky *create_Cholesky(double *A, int size);

/**
 * @brief Performs Cholesky decomposition on a symmetric positive-definite matrix.
 *
 * This function decomposes a symmetric positive-definite matrix A into a lower triangular matrix L
 * such that \( A = L \cdot L^T \). The resulting matrices L and Lᵀ are stored in the Cholesky structure.
 *
 * The function checks if the input matrix is symmetric and positive definite before performing the decomposition.
 *
 * @param A Pointer to the input square matrix (size: size x size).
 * @param size Dimension of the square matrix A (must be positive).
 *
 * @return Pointer to the Cholesky structure containing matrices A, L, and Lᵀ on success,
 *         or NULL on failure due to invalid dimensions, non-symmetric matrix, non-positive definite matrix,
 *         or memory allocation errors.
 */

Cholesky *Cholesky_decomposition(double *A, int size);

/**
 * @brief Frees all memory associated with a Cholesky decomposition structure.
 *
 * This function releases the memory allocated for the matrices A, L, and Lᵀ,
 * as well as the Cholesky structure itself.
 *
 * @param Cholesky_decomposition Pointer to the Cholesky structure to free.
 */

void free_Cholesky(Cholesky *Cholesky_decomposition);
    
/* generate_matrix.c */

/**
 * @brief Generates a matrix with random double values between 0 and DOUBLE.
 *
 * This function creates a matrix of size rows x columns, where each element
 * is a random double value in the range [0, DOUBLE].
 *
 * @param rows Number of rows in the matrix (must be positive).
 * @param columns Number of columns in the matrix (must be positive).
 *
 * @return Pointer to the generated matrix on success, or NULL on failure due to invalid dimensions
 *         or memory allocation errors.
 *
 */

double *generate_matrix_double (int rows, int columns);

/**
 * @brief Generates an identity matrix of given dimension.
 *
 * This function creates a square identity matrix of size dimension x dimension,
 * where all diagonal elements are 1.0 and all off-diagonal elements are 0.0.
 *
 * @param dimension Dimension of the identity matrix (must be positive).
 *
 * @return Pointer to the generated identity matrix on success, or NULL on failure due to invalid dimension
 *         or memory allocation errors.
 */

double *generate_identity_matrix (int dimension);

/* sequential_matrix_product.c */

/**
 * @brief Computes the product of two matrices P and Q sequentially.
 *
 * This function calculates the matrix product C = P * Q, where P is of size
 * P_rows x P_columns and Q is of size Q_rows x Q_columns. The computation is
 * performed sequentially without parallelization.
 *
 * @param P Pointer to the first input matrix (size: P_rows x P_columns).
 * @param P_rows Number of rows in matrix P (must be positive).
 * @param P_columns Number of columns in matrix P (must equal Q_rows).
 * @param Q Pointer to the second input matrix (size: Q_rows x Q_columns).
 * @param Q_rows Number of rows in matrix Q (must equal P_columns).
 * @param Q_columns Number of columns in matrix Q (must be positive).
 *
 * @return Pointer to the resulting matrix (size: P_rows x Q_columns) on success,
 *         or NULL on failure due to invalid dimensions, null pointers, or memory allocation errors.
 */

double *sequential_matrix_product(double *P, int P_rows, int P_columns, double *Q, int Q_rows, int Q_columns);

/* sequential_vector_matrix_product.c */

/**
 * @brief Computes the product of a matrix A and a vector X sequentially.
 *
 * This function calculates the vector result = A * X, where A is a matrix
 * and X is a vector. The computation is performed sequentially without parallelization.
 *
 * @param A Pointer to the input matrix (size: A_rows x A_columns).
 * @param A_rows Number of rows in the matrix A (must be positive).
 * @param A_columns Number of columns in the matrix A (must equal dimension).
 * @param X Pointer to the input vector (size: dimension).
 * @param dimension Dimension of the vector X (must equal A_columns).
 *
 * @return Pointer to the resulting vector (size: A_rows) on success,
 *         or NULL on failure due to invalid dimensions, null pointers, or memory allocation errors.
 */

double *sequential_vector_matrix_product(double *A, int A_rows, int A_columns, double *X, int dimension); 

/* parallel_vector_matrix_product.c */

/**
 * @brief Computes the product of a matrix A and a vector X in parallel using OpenMP.
 *
 * This function calculates the vector result = A * X, where A is a matrix
 * and X is a vector. The computation is parallelized using OpenMP to improve performance.
 *
 * @param A Pointer to the input matrix (size: A_rows x A_columns).
 * @param A_rows Number of rows in the matrix A (must be positive).
 * @param A_columns Number of columns in the matrix A (must equal dimension).
 * @param X Pointer to the input vector (size: dimension).
 * @param dimension Dimension of the vector X (must equal A_columns).
 *
 * @return Pointer to the resulting vector (size: A_rows) on success,
 *         or NULL on failure due to invalid dimensions, null pointers, or memory allocation errors.
 */

double *parallel_vector_matrix_product(double *A, int A_rows, int A_columns, double *X, int dimension);

/* parallel_matrix_product.c */

/**
 * @brief Computes the product of two matrices P and Q in parallel using OpenMP.
 *
 * This function calculates the matrix product C = P * Q using parallelization
 * to improve performance. Each element of the resulting matrix is computed
 * by summing the products of corresponding elements from rows of P and columns of Q.
 *
 * @param P Pointer to the first input matrix (size: P_rows x P_columns).
 * @param P_rows Number of rows in matrix P (must be positive).
 * @param P_columns Number of columns in matrix P (must equal Q_rows).
 * @param Q Pointer to the second input matrix (size: Q_rows x Q_columns).
 * @param Q_rows Number of rows in matrix Q (must equal P_columns).
 * @param Q_columns Number of columns in matrix Q (must be positive).
 *
 * @return Pointer to the resulting matrix (size: P_rows x Q_columns) on success,
 *         or NULL on failure due to invalid dimensions, null pointers, or memory allocation errors.
 */

double *parallel_matrix_product(double *P, int P_rows, int P_columns, double *Q, int Q_rows, int Q_columns);

/* vector_operations.c */

/**
 * @brief Computes the element-wise addition of two vectors X and Y.
 *
 * This function calculates the vector result = X + Y, where X and Y are vectors
 * of the same dimension.
 *
 * @param X Pointer to the first input vector (size: dimension).
 * @param Y Pointer to the second input vector (size: dimension).
 * @param dimension Dimension of the vectors (must be positive).
 *
 * @return Pointer to the resulting vector on success, or NULL on failure due to invalid dimension,
 *         null pointers, or memory allocation errors.
 */

double *vectors_addition(double *X, double *Y, int dimension);

/**
 * @brief Computes the element-wise product of two vectors X and Y.
 *
 * This function calculates the vector result = X .* Y, where .* denotes
 * element-wise multiplication.
 *
 * @param X Pointer to the first input vector (size: dimension).
 * @param Y Pointer to the second input vector (size: dimension).
 * @param dimension Dimension of the vectors (must be positive).
 *
 * @return Pointer to the resulting vector on success, or NULL on failure due to invalid dimension,
 *         null pointers, or memory allocation errors.
 */

double *scalar_product(double *X, double *Y, int dimension);

/**
 * @brief Computes the cross product of two 3-dimensional vectors X and Y.
 *
 * This function calculates the vector result = X x Y, where x denotes
 * the cross product. It assumes that both input vectors are 3-dimensional.
 *
 * @param X Pointer to the first input vector (size: 3).
 * @param Y Pointer to the second input vector (size: 3).
 *
 * @return Pointer to the resulting vector (size: 3) on success, or NULL on failure due to null pointers
 *         or memory allocation errors.
 */

double *vector_product(double *X, double *Y);

/**
 * @brief Computes the Euclidean norm (length) of a given vector X.
 *
 * This function calculates the norm ||X|| = sqrt(X_1^2 + X_2^2 + ... + X_n^2),
 * where n is the dimension of the vector.
 *
 * @param X Pointer to the input vector (size: dimension).
 * @param dimension Dimension of the vector (must be positive).
 *
 * @return The Euclidean norm of the vector on success, or -1.0 on failure due to invalid dimension
 *         or null pointers.
 */

double vector_norm(double *X, int dimension);

/* matrix_operations.c */

/**
 * @brief Adds two matrices A and B of the same size (rows x columns).
 *
 * This function computes the element-wise sum of two matrices A and B
 * and stores the result in a new matrix.
 *
 * @param A Pointer to the first input matrix (size: rows x columns).
 * @param B Pointer to the second input matrix (size: rows x columns).
 * @param rows Number of rows in the matrices (must be positive).
 * @param columns Number of columns in the matrices (must be positive).
 *
 * @return Pointer to the resulting matrix on success, or NULL on failure due to invalid dimensions,
 *         null pointers, or memory allocation errors.
 */

double *matrices_addition(double *A, double *B, int rows, int columns);

/**
 * @brief Multiplies a matrix A by a scalar value.
 *
 * This function multiplies each element of a matrix A by a given scalar value
 * and stores the result in a new matrix.
 *
 * @param A Pointer to the input matrix (size: rows x columns).
 * @param rows Number of rows in the matrix (must be positive).
 * @param columns Number of columns in the matrix (must be positive).
 * @param scalar Scalar value to multiply each element by.
 *
 * @return Pointer to the resulting matrix on success, or NULL on failure due to invalid dimensions,
 *         null pointers, or memory allocation errors.
 */

double *matrix_scalar_multiplication(double *A, int rows, int columns, int scalar);

/**
 * @brief Computes the transpose of a given matrix A.
 *
 * This function creates a new matrix that is the transpose of the input matrix A.
 *
 * @param A Pointer to the input matrix (size: rows x columns).
 * @param rows Number of rows in the input matrix (must be positive).
 * @param columns Number of columns in the input matrix (must be positive).
 *
 * @return Pointer to the transposed matrix on success, or NULL on failure due to invalid dimensions,
 *         null pointers, or memory allocation errors.
 */

double *matrix_transpose(double *A, int rows, int columns);

/**
 * @brief Computes the trace of a square matrix A.
 *
 * The trace is defined as the sum of the diagonal elements of a square matrix.
 *
 * @param A Pointer to the input square matrix (size: dimension x dimension).
 * @param dimension Dimension of the square matrix (must be positive).
 *
 * @return The trace of the matrix on success, or -1.0 on failure due to invalid dimensions
 *         or null pointers.
 */

double matrix_trace(double *A, int dimension);

/**
 * @brief Computes the 1-norm of a given matrix A.
 *
 * The 1-norm is defined as the maximum absolute column sum of a matrix.
 *
 * @param A Pointer to the input matrix (size: rows x columns).
 * @param rows Number of rows in the input matrix (must be positive).
 * @param columns Number of columns in the input matrix (must be positive).
 *
 * @return The 1-norm of the matrix on success, or -1.0 on failure due to invalid dimensions,
 *         null pointers, or memory allocation errors during transposition.
 */

double matrix_norm(double *A, int rows, int columns);

/**
 * @brief Computes the Frobenius norm of a given matrix A.
 *
 * The Frobenius norm is defined as the square root of the sum of squares
 * of all elements in a matrix.
 *
 * @param A Pointer to the input matrix (size: rows x columns).
 * @param rows Number of rows in the input matrix (must be positive).
 * @param columns Number of columns in the input matrix (must be positive).
 *
 * @return The Frobenius norm of the matrix on success, or -1.0 on failure due to invalid dimensions
 *         or null pointers.
 */

double frobenius_norm(double *A, int rows, int columns);

/**
 * @brief Computes the determinant of a square matrix A using LU decomposition.
 *
 * This function calculates the determinant by performing LU decomposition
 * and multiplying diagonal elements from L and U matrices. It checks for singularity during computation.
 *
 * @param A Pointer to the input square matrix (size: rows x rows).
 * @param rows Number of rows in the square matrix (must equal columns and be positive).
 * @param columns Number of columns in the square matrix (must equal rows and be positive).
 *
 * @return The determinant of the matrix on success, or -1.0 on failure due to invalid dimensions,
 *         null pointers, singularity detection, or LU decomposition failure.
 */

double matrix_determinant(double *A, int rows, int columns);

/**
 * @brief Computes eigenvalues of a square matrix using QR decomposition with iterative refinement.
 *
 * This function uses QR decomposition iteratively until convergence is achieved within
 * a specified tolerance or until reaching a maximum number of iterations. Eigenvalues are extracted from
 * diagonal elements after convergence.
 *
 * @param A Pointer to the input square matrix (size: rows x columns).
 * @param rows Number of rows in the input matrix (must be positive).
 * @param columns Number of columns in the input matrix (must be positive).
 * @param max_iter Maximum number of iterations for convergence.
 * @param tol Tolerance for convergence check (must be positive).
 *
 * @return Pointer to an array containing eigenvalues on success, or NULL on failure due to invalid dimensions,
           null pointers, memory allocation errors, or lack of convergence within max_iter iterations.
 */

double *matrix_eigenvalues(double *A, int rows, int columns, int max_iter, double tol);

/**
 * @brief Solves a lower triangular system Lc = b using forward substitution.
 *
 * This function computes the solution vector c such that Lc = b for a lower triangular matrix L.
 *
 * @param L Pointer to the lower triangular matrix (size: d x d).
 * @param d Dimension of the square matrix L and vector b (must be positive).
 * @param b Pointer to the right-hand side vector b (size: d).
 *
 * @return Pointer to the solution vector c on success, or NULL on failure due to invalid dimensions,
 *         null pointers, singular matrix detection, or memory allocation errors.
 */

double *forward_substitution(double *L, int d, double *b);

/**
 * @brief Solves an upper triangular system Ux = c using backward substitution.
 *
 * This function computes the solution vector x such that Ux = c for an upper triangular matrix U.
 *
 * @param U Pointer to the upper triangular matrix (size: d x d).
 * @param d Dimension of the square matrix U and vector c (must be positive).
 * @param c Pointer to the right-hand side vector c (size: d).
 *
 * @return Pointer to the solution vector x on success, or NULL on failure due to invalid dimensions,
 *         null pointers, singular matrix detection, or memory allocation errors.
 */

double *backward_substitution(double *U, int d, double *c);

/**
 * @brief Solves a linear system Ax = b using LU decomposition.
 *
 * This function computes the solution vector x for a square matrix A and right-hand side vector b
 * by performing LU decomposition and solving Ly = b followed by Ux = y.
 *
 * @param A Pointer to the square matrix A (size: d x d).
 * @param b Pointer to the right-hand side vector b (size: d).
 * @param d Dimension of the square matrix A and vector b (must be positive).
 *
 * @return Pointer to the solution vector x on success, or NULL on failure due to invalid dimensions,
 *         null pointers, LU decomposition failure, or memory allocation errors.
 */

double *solve_LU_system(double *A, double *b, int d);

/**
 * @brief Computes the inverse of a square matrix A using LU decomposition.
 *
 * This function calculates the inverse of a square matrix A by solving multiple systems
 * Ax_i = e_i, where e_i are columns of the identity matrix. The result is stored in a new matrix.
 *
 * @param A Pointer to the square matrix A (size: d x d).
 * @param d Dimension of the square matrix A (must be positive).
 *
 * @return Pointer to the inverse matrix on success, or NULL on failure due to invalid dimensions,
 *         null pointers, LU decomposition failure, or memory allocation errors.
 */

double *matrix_inverse(double *A, int n);

/**
 * @brief Checks if a square matrix H has converged based on a given tolerance.
 *
 * This function verifies whether all elements below the main diagonal are smaller than
 * a specified tolerance value. It assumes H is a square matrix.
 *
 * @param H Pointer to the square matrix (size: n x n).
 * @param n Dimension of the square matrix (must be positive).
 * @param tol Tolerance for convergence check (must be positive).
 *
 * @return 1 if converged, 0 if not converged, or -1 on failure due to invalid dimensions,
 *         null pointers, or invalid tolerance values.
 */

int has_converged(double *H, int n, double tol);

double *matrix_eigenvectors(double *A, int size, double *eigenvalues, int max_iter, double tol);

/* LU_decomposition.c */

/**
 * @brief Allocates and initializes an LU decomposition structure.
 *
 * This function creates an LU decomposition structure for a given matrix A
 * and initializes its components (L, U, and A).
 *
 * @param A Pointer to the input matrix (size: rows x columns).
 * @param rows Number of rows in the matrix (must be positive).
 * @param columns Number of columns in the matrix (must be positive).
 *
 * @return Pointer to the LU structure containing L, U, and A matrices on success,
 *         or NULL on failure due to invalid dimensions or memory allocation errors.
 */

LU *create_LU(double *A, int rows, int columns);

/**
 * @brief Performs LU decomposition on a given matrix using Doolittle's algorithm.
 *
 * This function decomposes a matrix A into a lower triangular matrix L
 * and an upper triangular matrix U such that A = L * U.
 *
 * @param A Pointer to the input matrix (size: rows x columns).
 * @param rows Number of rows in the matrix (must be positive).
 * @param columns Number of columns in the matrix (must be positive).
 *
 * @return Pointer to the LU structure containing L and U matrices on success,
 *         or NULL on failure due to invalid dimensions or memory allocation errors.
 */

LU *LU_decomposition(double *A, int rows, int columns);

/**
 * @brief Performs parallelized LU decomposition on a given matrix using OpenMP.
 *
 * This function decomposes a matrix A into a lower triangular matrix L
 * and an upper triangular matrix U such that A = L * U. The computation is
 * parallelized using OpenMP for improved performance on large matrices.
 *
 * @param A Pointer to the input matrix (size: rows x columns).
 * @param rows Number of rows in the matrix (must be positive).
 * @param columns Number of columns in the matrix (must be positive).
 *
 * @return Pointer to the LU structure containing L and U matrices on success,
 *         or NULL on failure due to invalid dimensions or memory allocation errors.
 */

LU *LU_decomposition_parallel(double *A, int rows, int columns);

/**
 * @brief Frees all memory associated with an LU decomposition structure.
 *
 * This function releases the memory allocated for the matrices L, U, and A,
 * as well as the LU structure itself.
 *
 * @param LU_decomposition Pointer to the LU structure to free.
 */

void LU_free(LU *LU_decomposition);

/* QR_decomposition.c */

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

QR *create_QR(double *A, int rows, int columns);

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

QR *QR_decomposition(double *A, int rows, int columns);

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

QR *QR_decomposition_parallel(double *A, int rows, int columns);

/**
 * @brief Frees all memory associated with a QR decomposition structure.
 *
 * This function releases the memory allocated for the matrices Q, R, and A,
 * as well as the QR structure itself.
 *
 * @param QR_decomposition Pointer to the QR structure to free.
 */

void QR_free(QR *QR_decomposition);

#endif 
