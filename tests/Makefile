FUNCTIONS_C = $(shell ls *.c)
CC = gcc
CFLAGS = -Wall -Werror
LIB = LinearAlgebraBasics.so

all : TEST_generate_matrix TEST_sequential_matrix_product
	./TEST_generate_matrix
	./TEST_sequential_matrix_product

TEST_generate_matrix : TEST_generate_matrix.c
	$(CC) $(CFLAGS) -o $@ $^ ./$(LIB)

TEST_sequential_matrix_product : TEST_sequential_matrix_product.c
	$(CC) $(CFLAGS) -o $@ $^ ./$(LIB)

clean :
	rm -f *.o *~
	rm -f LinearAlgebraBasics.so LinearAlgebraBasics.h
	rm -f TEST_generate_matrix TEST_sequential_matrix_product
