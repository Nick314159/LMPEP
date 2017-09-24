#!/bin/bash
#compile test program to executable file
gfortran -O3 test.f90 dslm1.f90 dslm2.f90 dslcorr1.f90 dzslcorr1.f90 dslcorr2.f90 dzslcorr2.f90 dsam.f90 dsacorr.f90 dzsacorr.f90 -llmpeptest -llmpep -o test.exe 

gfortran -O3 poly_test.f90 dslm1.f90 dslm2.f90 dslcorr1.f90 dzslcorr1.f90 dslcorr2.f90 dzslcorr2.f90 dsam.f90 dsacorr.f90 dzsacorr.f90 -llmpeptest -llmpep -o poly_test.exe 


