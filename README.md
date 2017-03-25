# LMPEP - Eigenvalue solver for general and specialized matrix polynomials #
This package is a collection of Fortran 90 subroutines for accurately 
and efficiently solving matrix eigenvalue problems using a variety of methods specfialized for distinct problem types.

## Authors ##
- [Thomas R. Cameron](http://thomasrcameron.com/), 
The College of Idaho, USA
- [Nikolas I. Steckley](http://www.nsteckley.com), 
Steckley & Associates


## Getting started ##
__LMPEP__ source code can be found in the _src_ directory. Numerical tests can be found in the _tests_ directory; please update the enviroment file in _tests/src/environment.f90_ with the path to the _tests_ directory to run them.

## Related articles ##
This software is based on the following articles:

 1. J. L. Aurentz, T. Mach, L. Robol, R. Vandebril, and D. S. Watkins, Fast and backward stable computation of the eigenvalues of matrix polynomials, Preprint on arXiv.org math, (2016).
 2. J. L. Aurentz, T. Mach, R. Vandebril, and D. S. Watkins, Fast and backward stable computation of roots of polynomials, SIAM J. Matrix Anal. Appl., 36 (2015), pp. 942–973.
 3. D. A. Bini, Numerical computation of polynomial zeros by means of Aberths method, Numer. Algorithms, 13 (1996), pp. 179–200.
 4. D. A. Bini, L. Gemignani, and F. Tisseur, The Ehrlich-Aberth method for the nonsymmetric tridiagonal eigenvalue problem, SIAM J. Matrix Anal. Appl., 27 (2005), pp. 153–175.
 5. D. A. Bini and V. Noferini, Solving polynomial eigenvalue problem by means of the Ehrlich-Aberth method, Linear Algebra Appl., 439 (2013), pp. 1130–1149.
 6. J. Gary, Hyman’s method applied to the general eigenvalue problem, Mathematics of Computation, 19 (1965), pp. 314–316.
 7. S. J. Hammerling, C. J. Munro, and T. Francoise, An algorithm for the complete solution of quadratic eigenvalue problem, Transactions on Mathematical Software, 39 (2013), p. 19.
 8. B. Parlett, Laguerre’s method applied to the matrix eigenvalue problem, Mathematics ofComputation, 18 (1964), pp. 464–485.
 9. B. Plestenjak, Numerical methods for the tridiagonal hyperbolic quadratic eigenvalue problem, SIAM J. Matrix Anal. Appl., 28 (2006), pp. 1157–1172.
 10. F. Tisseur, Backward error and condition of polynomial eigenvalue problem, Linear Algebra Appl., 309 (2000), pp. 339–361.
