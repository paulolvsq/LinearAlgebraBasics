FUNCTIONS_C = $(shell ls *.c)
FUNCTIONS_H = $(shell ls *.h)
FUNCTIONS_O = $(shell ls *.o)
CC = gcc
CFLAGS = -Wall -Werror

all : LinearAlgebraBasics.so LinearAlgebraBasics.h
	cp $^ ..
	cp $^ ../tests

LinearAlgebraBasics.so : generate_matrix.o sequential_matrix_product.o sequential_vector_matrix_product.o
	$(CC) -shared $(FUNCTIONS_O) -o $@

generate_matrix.o :
	$(CC) -c -o $@ generate_matrix.c

sequential_matrix_product.o :
	$(CC) -c -o $@ sequential_matrix_product.c

sequential_vector_matrix_product.o :
	$(CC) -c -o $@ sequential_vector_matrix_product.c

clean :
	rm -f *.o *~ 
	rm -f LinearAlgebraBasics.so
