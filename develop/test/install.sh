#!/bin/bash

#compile scalar subroutines to object files
gfortran -c -O3 balance.f90 DCFD.f90 DGR.f90 DRANDPOLYJT.f90 RESCHECK.f90 DAMVW.f90 DCFT.f90 DGTO2.f90 DCB.f90 DFCC.f90 DMQF.f90 DCDB.f90 DFGR.f90 DNORMALPOLY.f90 pzeros.f90

#add object files to library and copy library to /usr/lib64
ar crv libsptest.a *.o
sudo cp libsptest.a /usr/lib64/

#compile test program to executable file
gfortran -O3 test.f90 -lsptest -llmpep -o test.exe

#clean up
rm *.o
rm *.a
rm *.mod
