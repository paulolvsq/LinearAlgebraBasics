OS := $(shell uname -s)

ifeq ($(OS), Darwin)
    GCC_VERSIONS := $(shell ls /opt/homebrew/bin/gcc-* 2>/dev/null | grep -E 'gcc-[0-9]+' || which gcc)
    GCC := $(shell echo $(GCC_VERSIONS) | tr ' ' '\n' | grep -E 'gcc-[0-9]+' | sort -V | tail -n 1)
    ifeq ($(GCC),)
        GCC := gcc 
    endif
    SDK_PATH := $(shell xcrun --show-sdk-path)
    CC = $(GCC)
    CFLAGS = -Wall -Werror -isysroot $(SDK_PATH) -fopenmp
    LDFLAGS = -lm -lgomp
else
    CC = gcc
    CFLAGS = -Wall -Werror -fopenmp
    LDFLAGS = -lm
endif

LIB = LinearAlgebraBasics.so

TESTS = TEST_generate_matrix TEST_sequential_matrix_product TEST_sequential_vector_matrix_product TEST_vector_operations TEST_matrix_operations TEST_LU_decomposition TEST_QR_decomposition TEST_matrix_inverse TEST_parallel_vector_matrix_product TEST_parallel_matrix_product TEST_Cholesky TEST_LDLT

all : $(TESTS)
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
	./TEST_Cholesky
	./TEST_LDLT

TEST_generate_matrix : TEST_generate_matrix.c
	$(CC) $(CFLAGS) -o $@ $< ./$(LIB) $(LDFLAGS)

TEST_sequential_matrix_product : TEST_sequential_matrix_product.c
	$(CC) $(CFLAGS) -o $@ $< ./$(LIB) $(LDFLAGS)

TEST_sequential_vector_matrix_product : TEST_sequential_vector_matrix_product.c
	$(CC) $(CFLAGS) -o $@ $< ./$(LIB) $(LDFLAGS)

TEST_vector_operations : TEST_vector_operations.c
	$(CC) $(CFLAGS) -o $@ $< ./$(LIB) $(LDFLAGS)

TEST_matrix_operations : TEST_matrix_operations.c
	$(CC) $(CFLAGS) -o $@ $< ./$(LIB) $(LDFLAGS)

TEST_LU_decomposition : TEST_LU_decomposition.c
	$(CC) $(CFLAGS) -o $@ $< ./$(LIB) $(LDFLAGS)

TEST_QR_decomposition : TEST_QR_decomposition.c
	$(CC) $(CFLAGS) -o $@ $< ./$(LIB) $(LDFLAGS)

TEST_matrix_inverse : TEST_matrix_inverse.c
	$(CC) $(CFLAGS) -o $@ $< ./$(LIB) $(LDFLAGS)

TEST_parallel_vector_matrix_product : TEST_parallel_vector_matrix_product.c
	$(CC) $(CFLAGS) -o $@ $< ./$(LIB) $(LDFLAGS)

TEST_parallel_matrix_product : TEST_parallel_matrix_product.c
	$(CC) $(CFLAGS) -o $@ $< ./$(LIB) $(LDFLAGS)

TEST_Cholesky : TEST_Cholesky.c
	$(CC) $(CFLAGS) -o $@ $< ./$(LIB) $(LDFLAGS)

TEST_LDLT : TEST_LDLT.c
	$(CC) $(CFLAGS) -o $@ $< ./$(LIB) $(LDFLAGS)

clean :
	rm -f *.o *~
	rm -f $(LIB) LinearAlgebraBasics.h
	rm -f $(TESTS)
