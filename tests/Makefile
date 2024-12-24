CC = gcc
CFLAGS = -Wall -Werror
LIB = LinearAlgebraBasics.so
MATH = -lm
OMPFLAGS = -fopenmp

all : TEST_generate_matrix TEST_sequential_matrix_product TEST_sequential_vector_matrix_product TEST_vector_operations TEST_matrix_operations TEST_LU_decomposition TEST_QR_decomposition TEST_matrix_inverse TEST_parallel_vector_matrix_product TEST_parallel_matrix_product
	./TEST_generate_matrix
	./TEST_sequential_matrix_product
	./TEST_sequential_vector_matrix_product
	./TEST_vector_operations
	./TEST_matrix_operations
	./TEST_LU_decomposition
	./TEST_QR_decomposition
	./TEST_matrix_inverse
	./TEST_parallel_vector_matrix_product
	./TEST_parallel_matrix_product

TEST_generate_matrix : TEST_generate_matrix.c
	$(CC) $(CFLAGS) -o $@ $^ ./$(LIB) $(OMPFLAGS)

TEST_sequential_matrix_product : TEST_sequential_matrix_product.c
	$(CC) $(CFLAGS) -o $@ $^ ./$(LIB) $(OMPFLAGS)

TEST_sequential_vector_matrix_product : TEST_sequential_vector_matrix_product.c
	$(CC) $(CFLAGS) -o $@ $^ ./$(LIB) $(OMPFLAGS)

TEST_vector_operations : TEST_vector_operations.c
	$(CC) $(CFLAGS) -o $@ $^ ./$(LIB) $(MATH) $(OMPFLAGS)

TEST_matrix_operations : TEST_matrix_operations.c
	$(CC) $(CFLAGS) -o $@ $^ ./$(LIB) $(OMPFLAGS)

TEST_LU_decomposition : TEST_LU_decomposition.c
	$(CC) $(CFLAGS) -o $@ $^ ./$(LIB) $(OMPFLAGS)

TEST_QR_decomposition : TEST_QR_decomposition.c
	$(CC) $(CFLAGS) -o $@ $^ ./$(LIB) $(OMPFLAGS)

TEST_matrix_inverse : TEST_matrix_inverse.c
	$(CC) $(CFLAGS) -o $@ $^ ./$(LIB) $(OMPFLAGS)

TEST_parallel_vector_matrix_product : TEST_parallel_vector_matrix_product.c
	$(CC) $(CFLAGS) -o $@ $^ ./$(LIB) $(OMPFLAGS)

TEST_parallel_matrix_product : TEST_parallel_matrix_product.c
	$(CC) $(CFLAGS) -o $@ $^ ./$(LIB) $(OMPFLAGS)

clean :
	rm -f *.o *~
	rm -f LinearAlgebraBasics.so LinearAlgebraBasics.h
	rm -f TEST_generate_matrix TEST_sequential_matrix_product TEST_sequential_vector_matrix_product
	rm -f TEST_vector_operations TEST_matrix_operations TEST_LU_decomposition TEST_QR_decomposition
	rm -f TEST_matrix_inverse TEST_parallel_vector_matrix_product TEST_parallel_matrix_product
