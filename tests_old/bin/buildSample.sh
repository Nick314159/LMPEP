#!/bin/bash
FLAGS='-O3'
while [ ! $# -eq 0 ]
do
    case "$1" in
        --help | -h)
            echo 'Builds sample tridiagonal test executable. Moves executable one directory up. Use -d  for debug flags. -h for this help message'
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
gfortran $FLAGS ../environment.f90 util.f90 dslmpep_subroutines_test.f90 dgtlmpep_subroutines_test.f90 QEP3D/qep3deacx.f90 QEP3D/qep3dea.f90 QEP3D/qep3dlag.f90 sample_tri_test_driver.f90 -lscalapack -llapack -lblas && rm *.mod
mv a.out ../bin/sampleTri.out
