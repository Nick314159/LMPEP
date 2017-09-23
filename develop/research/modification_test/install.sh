#!/bin/bash

#compile all files in other_src
#gfortran -c -O3 other_src/*

#add object files to library and copy library to /usr/lib64
ar crv liblmpeptest.a *.o
#sudo cp liblmpeptest.a /usr/lib64/
sudo cp liblmpeptest.a /usr/lib/

#compile test program to executable file
gfortran -O3 test.f90 -llmpeptest -llmpep -o test.exe

#clean up
rm *.mod
rm *.o
rm *.a
