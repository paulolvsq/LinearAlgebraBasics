CC = gcc
CFLAGS = -Wall -Werror
LIB = LinearAlgebraBasics.so
MATH = -lm
OMPFLAGS = -fopenmp

all : PERF_LU_decomposition
	./PERF_LU_decomposition

PERF_LU_decomposition : PERF_LU_decomposition.c
	$(CC) $(CFLAGS) -o $@ $^ ./$(LIB) $(OMPFLAGS)

clean :
	rm -f *.o *~
	rm -f LinearAlgebraBasics.so LinearAlgebraBasics.h
	rm -f PERF_LU_decomposition
