all : execute

execute :
	cd functions && make
	cd tests && make

clean :
	rm -f *.o *~ *.so LinearAlgebraBasics.so LinearAlgebraBasics.h
	cd functions && make clean
	cd tests && make clean

.PHONY : all execute clean
