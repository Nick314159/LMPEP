#!/bin/bash

gfortran -c -O3 *.f90

ar crv libmodtest.a *.o
sudo cp libmodtest.a /usr/local/lib

rm *.o
rm *.a
