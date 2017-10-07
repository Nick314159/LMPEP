#!/bin/bash

#compile all files in other_src
gfortran -c -O3 other_src/*

#add object files to library and copy library to /usr/local/lib
ar crv liblmpeptest.a *.o
sudo cp liblmpeptest.a /usr/local/lib

#compile test program to executable file
gfortran -O3 test_dslm.f90 -L/usr/local/lib -llmpeptest -llmpep -o test_dslm.exe
gfortran -O3 dslm_comparisons.f90 -L/usr/local/lib -llmpeptest -llmpep -o dslm_comparisons.exe
gfortran -O3 poly_test.f90 -L/usr/local/lib -llmpeptest -llmpep -o poly_test.exe

#clean up
rm *.mod
rm *.o
rm *.a
