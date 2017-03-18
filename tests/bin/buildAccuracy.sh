#!/bin/bash
FLAGS='-O3'
while [ ! $# -eq 0 ]
do
    case "$1" in
        --help | -h)
            echo 'Builds spoly test executable. Moves executable one directory up. Use -d  for debug flags. -h for this help message'
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
gfortran $FLAGS ../environment.f90 util.f90 dslmpep_subroutines_test.f90 dgelmpep_subroutines_test.f90 quadratic-eigensolver/dg3evx/izmaxa.f90 quadratic-eigensolver/dg3evx/dselct.f90 quadratic-eigensolver/dg3evx/dlasge.f90 quadratic-eigensolver/dg3evx/dlaqp3.f90 quadratic-eigensolver/dg3evx/dlanab.f90 quadratic-eigensolver/dg3evx/dlag3v.f90 quadratic-eigensolver/dg3evx/dlag3r.f90 quadratic-eigensolver/dg3evx/dlag3b.f90 quadratic-eigensolver/dg3evx/dlag3c.f90 quadratic-eigensolver/dg3evx/dg3evx.f90 accuracy_test_driver.f90 -lscalapack -llapack -lblas && rm *.mod
mv a.out ../bin/accuracy.out
