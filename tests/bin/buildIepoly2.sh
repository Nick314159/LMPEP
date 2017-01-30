#!/bin/bash
FLAGS='-O3'
while [ ! $# -eq 0 ]
do
    case "$1" in
        --help | -h)
            echo 'Builds iepoly2 test executable. Moves executable one directory up. Use -d  for debug flags. -h for this help message'
            exit
            ;;
        --debug | -d)
            echo 'Building with debug...' 
            FLAGS=$FLAGS' -g -Wall -Wextra'
            ;;
    esac
    shift
done
cd ../src
gfortran $FLAGS environment.f90 util.f90 ../../src/dslmpep_subroutines.f90 ../../src/dgelmpep_subroutines.f90 ie_test_driver2.f90 -lscalapack-openmpi -llapack -lblas && rm *.mod
mv a.out ../bin/iepoly2.out
