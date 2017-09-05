#!/bin/bash

#compile all files in lmpep_src
gfortran -c -O3 lmpep_src/*

#add object files to library and copy library to /usr/lib64
ar crv liblmpep.a *.o
sudo cp liblmpep.a /usr/lib64/

#clean up
rm *.o
rm *.a

#Generate documentation
#doxygen Doxygen/config
