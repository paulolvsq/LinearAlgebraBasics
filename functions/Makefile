OS := $(shell uname -s)

ifeq ($(OS), Darwin)
    GCC_VERSIONS := $(shell ls /opt/homebrew/bin/gcc-* 2>/dev/null | grep -E 'gcc-[0-9]+' || which gcc)
    GCC := $(shell echo $(GCC_VERSIONS) | tr ' ' '\n' | grep -E 'gcc-[0-9]+' | sort -V | tail -n 1)
    ifeq ($(GCC),)
        GCC := gcc  
    endif
    SDK_PATH := $(shell xcrun --show-sdk-path)
    CC = $(GCC)
    CFLAGS = -Wall -Werror -isysroot $(SDK_PATH) -fopenmp -lm
    OMPFLAGS = -lgomp
else
    CC = gcc
    CFLAGS = -Wall -Werror -fopenmp -fPIC -lm
    OMPFLAGS = 
endif

LIB = LinearAlgebraBasics.so

OBJECTS = generate_matrix.o sequential_matrix_product.o sequential_vector_matrix_product.o vector_operations.o matrix_operations.o LU_decomposition.o QR_decomposition.o parallel_vector_matrix_product.o parallel_matrix_product.o Cholesky_decomposition.o LDLT_decomposition.o

all : $(LIB) LinearAlgebraBasics.h
	cp $^ ..
	cp $^ ../tests
	cp $^ ../performances

$(LIB) : $(OBJECTS)
	$(CC) $(CFLAGS) -shared $^ -o $@

generate_matrix.o : generate_matrix.c
	$(CC) $(CFLAGS) -c -o $@ $<

sequential_matrix_product.o : sequential_matrix_product.c
	$(CC) $(CFLAGS) -c -o $@ $<

sequential_vector_matrix_product.o : sequential_vector_matrix_product.c
	$(CC) $(CFLAGS) -c -o $@ $<

vector_operations.o : vector_operations.c
	$(CC) $(CFLAGS) -c -o $@ $<

matrix_operations.o : matrix_operations.c
	$(CC) $(CFLAGS) -c -o $@ $<

Cholesky_decomposition.o : Cholesky_decomposition.c
	$(CC) $(CFLAGS) -c -o $@ $<

LDLT_decomposition.o : LDLT_decomposition.c
	$(CC) $(CFLAGS) -c -o $@ $<

parallel_matrix_product.o : parallel_matrix_product.c
	$(CC) $(CFLAGS) $(OMPFLAGS) -c -o $@ $<

parallel_vector_matrix_product.o : parallel_vector_matrix_product.c
	$(CC) $(CFLAGS) $(OMPFLAGS) -c -o $@ $<

LU_decomposition.o : LU_decomposition.c
	$(CC) $(CFLAGS) $(OMPFLAGS) -c -o $@ $<

QR_decomposition.o : QR_decomposition.c
	$(CC) $(CFLAGS) $(OMPFLAGS) -c -o $@ $<

clean :
	rm -f *.o *~ 
	rm -f $(LIB)
