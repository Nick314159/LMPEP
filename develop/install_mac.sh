#!/bin/bash

#compile random generator subroutines
gfortran -c -O3 init_random_seed.f90 daruv.f90

#compile Bini subroutines
gfortran -c -O3 cmerge.f90 cnvex.f90 ctest.f90 left.f90 right.f90

#compile scalar LMPEP subroutines
gfortran -c -O3 dzmod.f90 dsstart.f90 dseval.f90 drevseval.f90 dzseval.f90 dzrevseval.f90 dslcorr.f90 dzslcorr.f90 dslm.f90

#add object files to library and copy library to /usr/lib64
ar crv liblmpep.a *.o
sudo cp liblmpep.a /usr/local/lib

#compile test program to executable file
gfortran -O3 test_sp.f90 -llmpep -o test_sp.exe

#clean up
rm *.o
rm *.a

#Generate documentation
#doxygen Doxygen/config