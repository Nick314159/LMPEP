#!/bin/bash
FLAGS='-O3'
while [ ! $# -eq 0 ]
do
    case "$1" in
        --help | -h)
            echo 'Builds Computational Complexity test executable. Moves executable one directory up. Use -d  for debug flags. -h for this help message'
            exit
            ;;
        --debug | -d)
            echo 'Building with debug...'
            FLAGS=$FLAGS' -g -Wall -Wextra'
            ;;
    esac
    shift
done

gfortran $FLAGS environment.f90 util.f90 dslmpep_subroutines.f90  dgelmpep_subroutines.f90 dgeeam_subroutines.f90 dhslmpep_subroutines.f90 cc_test_driver.f90 -lscalapack-openmpi -llapack -lblas && rm *.mod
mv a.out ../bin/complexity.out
