FUNCTIONS_C = $(shell ls *.c)
FUNCTIONS_H = $(shell ls *.h)
FUNCTIONS_O = $(shell ls *.o)
CC = gcc
CFLAGS = -Wall -Werror
MATH = -lm

all : LinearAlgebraBasics.so LinearAlgebraBasics.h
	cp $^ ..
	cp $^ ../tests

LinearAlgebraBasics.so : generate_matrix.o sequential_matrix_product.o sequential_vector_matrix_product.o vector_operations.o matrix_operations.o LU_decomposition.o QR_decomposition.o
	$(CC) -shared $(FUNCTIONS_O) -o $@ $(MATH)

generate_matrix.o :
	$(CC) -c -o $@ generate_matrix.c

sequential_matrix_product.o :
	$(CC) -c -o $@ sequential_matrix_product.c

sequential_vector_matrix_product.o :
	$(CC) -c -o $@ sequential_vector_matrix_product.c

vector_operations.o : 
	$(CC) -c -o $@ vector_operations.c $(MATH)

matrix_operations.o :
	$(CC) -c -o $@ matrix_operations.c $(MATH)

LU_decomposition.o :
	$(CC) -c -o $@ LU_decomposition.c

QR_decomposition.o :
	$(CC) -c -o $@ QR_decomposition.c $(MATH)

clean :
	rm -f *.o *~ 
	rm -f LinearAlgebraBasics.so
