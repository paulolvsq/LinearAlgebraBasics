FUNCTIONS_C = $(shell ls *.c)
CC = gcc
CFLAGS = -Wall -Werror
LIB = LinearAlgebraBasics.so
MATH = -lm

all : TEST_generate_matrix TEST_sequential_matrix_product TEST_sequential_vector_matrix_product TEST_vector_operations TEST_matrix_operations TEST_LU_decomposition
	./TEST_generate_matrix
	./TEST_sequential_matrix_product
	./TEST_sequential_vector_matrix_product
	./TEST_vector_operations
	./TEST_matrix_operations
	./TEST_LU_decomposition

TEST_generate_matrix : TEST_generate_matrix.c
	$(CC) $(CFLAGS) -o $@ $^ ./$(LIB)

TEST_sequential_matrix_product : TEST_sequential_matrix_product.c
	$(CC) $(CFLAGS) -o $@ $^ ./$(LIB)

TEST_sequential_vector_matrix_product : TEST_sequential_vector_matrix_product.c
	$(CC) $(CFLAGS) -o $@ $^ ./$(LIB)

TEST_vector_operations : TEST_vector_operations.c
	$(CC) $(CFLAGS) -o $@ $^ ./$(LIB) $(MATH)

TEST_matrix_operations : TEST_matrix_operations.c
	$(CC) $(CFLAGS) -o $@ $^ ./$(LIB) 

TEST_LU_decomposition : TEST_LU_decomposition.c
	$(CC) $(CFLAGS) -o $@ $^ ./$(LIB) 

clean :
	rm -f *.o *~
	rm -f LinearAlgebraBasics.so LinearAlgebraBasics.h
	rm -f TEST_generate_matrix TEST_sequential_matrix_product TEST_sequential_vector_matrix_product
	rm -f TEST_vector_operations TEST_matrix_operations
