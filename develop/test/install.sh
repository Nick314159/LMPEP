#!/bin/bash

#compile scalar subroutines to object files
gfortran -c -O3 balance.f90 DCFD.f90 DGR.f90 DRANDPOLYJT.f90 RESCHECK.f90 DAMVW.f90 DCFT.f90 DGTO2.f90 init_random_seed.f90 DCB.f90 DFCC.f90 DMQF.f90 DCDB.f90 DFGR.f90 DNORMALPOLY.f90 pzeros.f90 pzeros_drive.f90

#add object files to library and copy library to /usr/lib64
ar crv libtest.a *.o
sudo cp libtest.a /usr/lib64/

#clean up
rm *.o
rm *.a
rm *.mod
