#!/bin/bash

gfortran -c -O3 src/*

ar crv libmodtest.a *.o
sudo cp libmodtest.a /usr/local/lib

gfortran -O3 time_and_error.f90 -L/usr/local/lib -lmodtest -llmpep -o time_and_error.exe
gfortran -O3 poly_test.f90 -L/usr/local/lib -lmodtest -llmpep -o poly_test.exe

rm *.o
rm *.a
