# Détection du système d'exploitation
OS := $(shell uname -s)

# Configuration spécifique à macOS
ifeq ($(OS), Darwin)
    GCC_VERSIONS := $(shell ls /opt/homebrew/bin/gcc-* 2>/dev/null | grep -E 'gcc-[0-9]+' || which gcc)
    GCC := $(shell echo $(GCC_VERSIONS) | tr ' ' '\n' | grep -E 'gcc-[0-9]+' | sort -V | tail -n 1)
    ifeq ($(GCC),)
        GCC := gcc  # Utiliser le gcc standard si aucun gcc spécifique n'est trouvé
    endif
    SDK_PATH := $(shell xcrun --show-sdk-path)
    CC = $(GCC)
    CFLAGS = -Wall -Werror -isysroot $(SDK_PATH) -fopenmp
    LDFLAGS = -lm -lgomp
else
    # Configuration standard pour Linux
    CC = gcc
    CFLAGS = -Wall -Werror -fopenmp
    LDFLAGS = -lm
endif

LIB = LinearAlgebraBasics.so

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
