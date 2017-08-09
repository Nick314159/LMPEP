#!/bin/bash

#compile scalar subroutines to object files
gfortran -c -O3 util.f90 dzmod.f90 dsstart.f90 dseval.f90 zseval.f90

#add object files to library and copy library to /usr/lib64
ar crv libtest.a *.o
sudo cp libtest.a /usr/lib64/

#compile test program to executable file
gfortran -O3 test.f90 -ltest -o test.exe

#clean up
rm *.mod
rm *.o
rm *.a
