#!/bin/bash
gfortran util.f90 dslmpep_subroutines.f90 dgelmpep_subroutines.f90 dgeeam_subroutines.f90 ie_test_driver2.f90 -lscalapack-openmpi -llapack -lblas && rm *.mod

