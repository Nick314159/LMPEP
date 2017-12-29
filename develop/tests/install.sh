#!/bin/bash

#compile all files in lmpep_src
gfortran -c -O3 src/*

#add object files to library and copy library to /usr/local/lib
ar crv liblmpeptest.a *.o
cp liblmpeptest.a /usr/local/lib

#compile testMethods.f90 file and create executable output testMethods
gfortran -O3 testMethods.f90 -L/usr/local/lib -llmpep -llmpeptest -o testMethods

#compile testConvergence.f90 file and create executable output testConvergence
gfortran -O3 testConvergence.f90 -L/usr/local/lib -llmpep -llmpeptest -o testConvergence

#compile testComparison.f90 file and create executable output testComparison
gfortran -O3 testComparison.f90 -L/usr/local/lib -llmpep -llmpeptest -lblas -o testComparison

#clean up
rm *.mod
rm *.o
rm *.a
