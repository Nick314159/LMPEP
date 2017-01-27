#!/bin/bash
FLAGS='-O3'
while [ ! $# -eq 0 ]
do
    case "$1" in
        --help | -h)
            echo 'Builds tridiagonal test executable. Moves executable one directory up. Use -d  for debug flags. -h for this help message'
            exit
            ;;
        --debug | -d)
            echo 'Building with debug...'
            FLAGS=$FLAGS' -g -Wall -Wextra'
            ;;
    esac
    shift
done

gfortran $FLAGS environment.f90 util.f90 dslmpep_subroutines.f90 dgtlmpep_subroutines.f90 eigen_v1.1/bgt.f90 tri_test_driver.f90 -lscalapack-openmpi -llapack -lblas && rm *.mod
mv a.out ../bin/tri.out
