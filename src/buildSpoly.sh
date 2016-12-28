#!/bin/bash
gfortran util.f90 dslmpep_subroutines.f90 pzeros.f90 spoly_test_driver.f90 -lscalapack-openmpi -llapack -lblas && rm *.mod

