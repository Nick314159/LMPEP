`quadratic-eigensolver` - MATLAB, Octave and Fortran codes for solving quadratic eigenvalue problems
====================================================================================================

About
-----
`quadratic-eigensolver' contains a MATLAB function, an Octave function and
Fortran routines for the numerical solution of quadratic eigenvalue problems
based on the algorithm in the paper:

S. Hammarling, C. J. Munro and F. Tisseur, "[An algorithm for the
complete solution of quadratic eigenvalue problems]
(http://eprints.ma.man.ac.uk/2061/)",
ACM Trans. Math. Software, 39(3):18:1-18:19, 2013.

Description of the main directories:

* MATLAB: contains the MATLAB function `quadeig`. By default, it uses the
  [NAG Toolbox for MATLAB] if available
  (http://www.nag.co.uk/numeric/MB/start.asp).

* dg3evx: this directory contains double-precision Fortran 90 files for real
  quadratic eigenvalue problems. The driver subroutine is `dg3evx`. The other
  Fortran files are auxiliary subroutines. The LAPACK and BLAS libraries are
  required to run `dg3evx`.

* zg3evx: this directory contains double-precision complex Fortran 90 files for
  complex quadratic eigenvalue problems. The main subroutine is `zg3evx`. The
  other Fortran files are auxiliary subroutines. The LAPACK and BLAS libraries
  are required to run `zg3evx`.

* Octave: contains the Octave function `quadeig_o`.

License
-------

See `license.txt` for licensing information.

