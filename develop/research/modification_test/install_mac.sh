#!/bin/bash

gfortran -c -O3 *.f90

ar crv libmodtest.a *.o
sudo cp libmodtest.a /usr/local/lib

gfortran -O3 test.f90 -L/usr/local/lib -lmodtest -llmpep -o test.exe

rm *.o
rm *.a
