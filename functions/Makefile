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
    CFLAGS = -Wall -Werror -isysroot $(SDK_PATH) -fopenmp -lm
    OMPFLAGS = -lgomp
else
    # Configuration standard pour Linux
    CC = gcc
    CFLAGS = -Wall -Werror -fopenmp -fPIC -lm
    OMPFLAGS = 
endif

LIB = LinearAlgebraBasics.so

all : LinearAlgebraBasics.so LinearAlgebraBasics.h
	cp $^ ..
	cp $^ ../tests
	cp $^ ../performances

LinearAlgebraBasics.so : generate_matrix.o sequential_matrix_product.o sequential_vector_matrix_product.o vector_operations.o matrix_operations.o LU_decomposition.o QR_decomposition.o parallel_vector_matrix_product.o parallel_matrix_product.o
	$(CC) $(CFLAGS) -shared $^ -o $@ $(MATH)

generate_matrix.o : generate_matrix.c
	$(CC) $(CFLAGS) -c -o $@ $^

sequential_matrix_product.o : sequential_matrix_product.c
	$(CC) $(CFLAGS) -c -o $@ $^ 

parallel_matrix_product.o : parallel_matrix_product.c
	$(CC) $(CFLAGS) $(OMPFLAGS) -c -o $@ $^

sequential_vector_matrix_product.o : sequential_vector_matrix_product.c
	$(CC) $(CFLAGS) -c -o $@ $^

parallel_vector_matrix_product.o : parallel_vector_matrix_product.c
	$(CC) $(CFLAGS) $(OMPFLAGS) -c -o $@ $^

vector_operations.o : vector_operations.c
	$(CC) $(CFLAGS) -c -o $@ $^ $(MATH)

matrix_operations.o : matrix_operations.c
	$(CC) $(CFLAGS) -c -o $@ $^ $(MATH)

LU_decomposition.o : LU_decomposition.c
	$(CC) $(CFLAGS) $(OMPFLAGS) -c -o $@ $^

QR_decomposition.o : QR_decomposition.c
	$(CC) $(CFLAGS) $(OMPFLAGS) -c -o $@ $^ $(MATH)

clean :
	rm -f *.o *~ 
	rm -f LinearAlgebraBasics.so
