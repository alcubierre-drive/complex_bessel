CC ?= gcc
LD ?= gcc
FF ?= gfortran

LIBS := -lgfortran -lm
CFLAGS := -fPIC -Wall -Wextra -pedantic -std=c11
FFLAGS := -fPIC

libbessel.so: bessel.o mod_zbes.o
	$(LD) $^ -o $@ -shared $(LIBS)

%.o: %.c
	$(CC) -c $^ -o $@ $(CFLAGS)
%.o: %.f90
	$(FF) -c $^ -o $@ $(FFLAGS)
