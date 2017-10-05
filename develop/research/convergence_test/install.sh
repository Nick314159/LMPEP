#!/bin/bash

gfortran -c -O3 src/* 

ar crv libmodtest.a *.o
sudo cp libmodtest.a /usr/local/lib

gfortran -O3 conv_test.f90 -L/usr/local/lib -lmodtest -llmpep -llapack -lblas  -o conv_test.exe

rm *.o
rm *.a
