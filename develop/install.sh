#!/bin/bash

#compile all files in lmpep_src
gfortran -c -O3 lmpep_src/*

#add object files to library and copy library to /usr/local/lib
ar crv liblmpep.a *.o
cp liblmpep.a /usr/local/lib

#clean up
rm *.o
rm *.a

#Generate documentation
doxygen Doxygen/config
