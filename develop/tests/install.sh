#!/bin/bash

gfortran -c -O3 src/*

ar crv liblmpeptest.a *.o
cp liblmpeptest.a /usr/local/lib

gfortran -O3 testMethods.f95 -L/usr/local/lib -llmpep -llmpeptest -o test

rm *.o
rm *.a
