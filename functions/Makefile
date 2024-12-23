CC = gcc
CFLAGS = -Wall -Werror -fPIC
MATH = -lm
OMPFLAGS = -fopenmp

all : LinearAlgebraBasics.so LinearAlgebraBasics.h
	cp $^ ..
	cp $^ ../tests
	cp $^ ../performances

LinearAlgebraBasics.so : generate_matrix.o sequential_matrix_product.o sequential_vector_matrix_product.o vector_operations.o matrix_operations.o LU_decomposition.o QR_decomposition.o 
	$(CC) $(CFLAGS) -shared $^ -o $@ $(MATH)

generate_matrix.o : generate_matrix.c
	$(CC) $(CFLAGS) -c -o $@ $^

sequential_matrix_product.o : sequential_matrix_product.c
	$(CC) $(CFLAGS) -c -o $@ $^ 

sequential_vector_matrix_product.o : sequential_vector_matrix_product.c
	$(CC) $(CFLAGS) -c -o $@ $^

vector_operations.o : vector_operations.c
	$(CC) $(CFLAGS) -c -o $@ $^ $(MATH)

matrix_operations.o : matrix_operations.c
	$(CC) $(CFLAGS) -c -o $@ $^ $(MATH)

LU_decomposition.o : LU_decomposition.c
	$(CC) $(CFLAGS) $(OMPFLAGS) -c -o $@ $^

QR_decomposition.o : QR_decomposition.c
	$(CC) $(CFLAGS) -c -o $@ $^ $(MATH)

clean :
	rm -f *.o *~ 
	rm -f LinearAlgebraBasics.so
