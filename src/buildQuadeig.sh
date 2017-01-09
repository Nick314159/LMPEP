#!/bin/bash
FLAGS='-O3'
while [ ! $# -eq 0 ]
do
    case "$1" in
        --help | -h)
            echo 'Builds iepoly1 test executable. Moves executable one directory up. Use -d  for debug flags. -h for this help message'
            exit
            ;;
        --debug | -d)
            echo 'Building with debug...'
            FLAGS=$FLAGS' -g -Wall -Wextra'
            ;;
    esac
    shift
done
cd quadratic-eigensolver/
gfortran $FLAGS ../environment.f90 ../util.f90 ../dslmpep_subroutines.f90 ../dgelmpep_subroutines.f90 ../dgeeam_subroutines.f90  dlag3r.f90 dlaqp3.f90 izmaxa.f90 dg3evx.f90 dlag3v.f90 dlasge.f90 dlag3c.f90 dlanab.f90 dselct.f90 quadeig_test_driver.f90  -lscalapack-openmpi -llapack -lblas && rm *.mod
mv a.out ../../bin/quadeig.out



