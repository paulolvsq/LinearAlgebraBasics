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
    CFLAGS = -Wall -Werror -g -isysroot $(SDK_PATH) -fopenmp
    LDFLAGS = -lm -lgomp
else
    # Configuration standard pour Linux
    CC = gcc
    CFLAGS = -Wall -Werror -g -fopenmp
    LDFLAGS = -lm
endif

LIB = LinearAlgebraBasics.so

all : PERF_LU_decomposition PERF_QR_decomposition PERF_vector_matrix_product PERF_matrix_product
	./PERF_LU_decomposition
	./PERF_QR_decomposition
	./PERF_vector_matrix_product
	./PERF_matrix_product

PERF_LU_decomposition : PERF_LU_decomposition.c
	$(CC) $(CFLAGS) -o $@ $^ ./$(LIB) $(LDFLAGS)

PERF_QR_decomposition : PERF_QR_decomposition.c
	$(CC) $(CFLAGS) -o $@ $^ ./$(LIB) $(LDFLAGS)

PERF_vector_matrix_product : PERF_vector_matrix_product.c 
	$(CC) $(CFLAGS) -o $@ $^ ./$(LIB) $(LDFLAGS)

PERF_matrix_product : PERF_matrix_product.c
	$(CC) $(CFLAGS) -o $@ $^ ./$(LIB) $(LDFLAGS)

clean :
	rm -f *.o *~
	rm -f LinearAlgebraBasics.so LinearAlgebraBasics.h
	rm -f PERF_LU_decomposition PERF_QR_decomposition PERF_vector_matrix_product PERF_matrix_product
