#!/bin/bash
#compile test program to executable file
gfortran -O3 test.f90 dslm1.f90 dslcorr1.f90 dzslcorr1.f90 dsam.f90 dsacorr.f90 dzsacorr.f90 -llmpeptest -llmpep -o test.exe

