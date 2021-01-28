.PHONY: all lib install clean

#F90=gfortran -fno-second-underscore -malign-double -ggdb -DCERNLIB_DOUBLE
F90=gfortran -fno-second-underscore -malign-double -O3 -DCERNLIB_DOUBLE
PATHLIB=./lib

all: lib

lib: 
	$(F90) -c *.F
	ar rcs libminuit.a *.o

install: all
	cp libminuit.a  $(PATHLIB)

clean: 
	rm -f *.o
	rm -f *.mod
	rm -f *.a
