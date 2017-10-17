#!/bin/bash

#compile all files in lmpep_src
gfortran -c -O3 src/*

#add object files to library and copy library to /usr/local/lib
ar crv liblmpeptest.a *.o
cp liblmpeptest.a /usr/local/lib

#compile test.f90 file and create executable output testMethods
gfortran -O3 testMethods.f90 -L/usr/local/lib -llmpep -llmpeptest -o testMethods

rm *.o
rm *.a
