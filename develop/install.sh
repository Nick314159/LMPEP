#!/bin/bash
gfortran -c -O3 *.f90
ar crv libtest.a *.o
sudo cp libtest.a /usr/lib64/
gfortran -O3 test.f90 -ltest -o test.sh
