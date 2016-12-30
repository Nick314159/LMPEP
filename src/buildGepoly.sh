#!/bin/bash
gfortran -g -Wall environment.f90 util.f90 dslmpep_subroutines.f90 dgelmpep_subroutines.f90 dgeeam_subroutines.f90 ge_test_driver.f90 -lscalapack-openmpi -llapack -lblas && rm *.mod

