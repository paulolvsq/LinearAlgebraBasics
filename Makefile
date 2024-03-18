MAKE = make

all : execute

execute :
	cd functions && $(MAKE)
	cd tests && $(MAKE)

clean :
	rm -f *.o *~ *.so LinearAlgebraBasics.so LinearAlgebraBasics.h
	cd functions && $(MAKE) clean
	cd tests && $(MAKE) clean

.PHONY : all execute clean
