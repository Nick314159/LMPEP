#gfortran  util.f90 dslmpep_subroutines.f90 pzeros.f90 test_driver.f90 -lscalapack-openmpi -llapack -lblas

#Define variables
f90comp = gfortran
switch = -O3
objects =   util.f90 dslmpep_subroutines.f90 pzeros.f90 test_driver.f90
dependencies = -lscalapack-openmpi -llapack -lblas
#Makefile
main.out: $(objects)
	$(f90comp) -o main.out $(switch) $(objects) $(dependecies)

util.mod: util.o util.f90
	$(f90comp) -c $(switch) util.f90 $(dependecies)

util.o: util.f90
	$(f90comp) -c $(switch) util.f90 $(dependecies)

dslmpep_subroutines.mod: dslmpep_subroutines.o dslmpep_subroutines.f90
	$(f90comp) -c $(switch) dslmpep_subroutines.f90 $(dependecies)

dslmpep_subroutines.o: dslmpep_subroutines.f90 dgelmpep_subroutines.mod
	$(f90comp) -c $(switch) dslmpep_subroutines.f90 $(dependecies)

poly_zero.mod: poly_zero.o poly_zero.f90
	$(f90comp) -c $(switch) poly_zero.f90 $(dependecies)

poly_zero.o: poly_zero.f90
	$(f90comp) -c $(switch) poly_zero.f90 $(dependecies)
	
test_driver.o: test_driver.f90
	$(f90comp) -c $(switch) test_driver.f90 $(dependecies)

clean:
	rm *.mod
	rm fort.1

