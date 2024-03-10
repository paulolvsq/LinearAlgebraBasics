FUNCTIONS_C = $(shell ls *.c)
FUNCTIONS_H = $(shell ls *.h)
FUNCTIONS_O = $(shell ls *.o)
CC = gcc
CFLAGS = -Wall -Werror

all : LinearAlgebraBasics.so LinearAlgebraBasics.h
	cp $^ ..
	cp $^ ../tests

LinearAlgebraBasics.so : generate_matrix.o
	$(CC) -shared $(FUNCTIONS_O) -o $@

generate_matrix.o :
	$(CC) -c -o $@ generate_matrix.c

clean :
	rm -f *.o *~
	rm -f LinearAlgebraBasics.so
