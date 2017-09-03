# LMPEP - Laguerre's Method Applied to the Polynomial Eigenvalue Problem ##
LMPEP is a Fortran library collection of solvers and their dependencies for the polynomial eigenvalue problem. There are several subroutines that make use of external subroutines from the BLAS and LAPACK libraries. 

## Authors ##
- [Thomas R. Cameron](http://thomasrcameron.com/), 
Davidson College, NC
- [Nikolas I. Steckley](http://www.nsteckley.com), 
Portland State University, OR


## In Progress ##
Currently we are updating the code in the develop section, in the hopes of developing more robust algorithms and better documentation. Also, we are working on the parallelization of these methods. 

## Related articles ##
This software is based on the following articles:

 1. J. L. Aurentz, T. Mach, L. Robol, R. Vandebril, and D. S. Watkins, Fast and backward stable computation of the eigenvalues of matrix polynomials, Preprint on arXiv.org math, (2016).
 2. J. L. Aurentz, T. Mach, R. Vandebril, and D. S. Watkins, Fast and backward stable computation of roots of polynomials, SIAM Journal on Matrix Analysis and Applications, 36(3) (2015), pp. 942–973.
 3. D. A. Bini, Numerical computation of polynomial zeros by means of Aberths method, Numerical Algorithms, 13(2) (1996), pp. 179–200.
 4. D. A. Bini, L. Gemignani, and F. Tisseur, The Ehrlich-Aberth method for the nonsymmetric tridiagonal eigenvalue problem, SIAM Journal on Matrix Analysis and Applications, 27(1) (2005), pp. 153–175.
 5. D. A. Bini and V. Noferini, Solving polynomial eigenvalue problem by means of the Ehrlich-Aberth method, Linear Algebra and its Applications, 439 (2013), pp. 1130–1149.
 6. J. Gary, Hyman’s method applied to the general eigenvalue problem, Mathematics of Computation, 19 (1965), pp. 314–316.
 7. S. J. Hammarling, C. J. Munro, and F. Tisseur, An algorithm for the complete solution of quadratic eigenvalue problem, Transactions on Mathematical Software, 39(3) (2013), p. 19.
 8. B. Parlett, Laguerre’s method applied to the matrix eigenvalue problem, Mathematics ofComputation, 18 (1964), pp. 464–485.
 9. B. Plestenjak, Numerical methods for the tridiagonal hyperbolic quadratic eigenvalue problem, SIAM Journal on Matrix Analysis and Applications, 28(4) (2006), pp. 1157–1172.
 10. F. Tisseur, Backward error and condition of polynomial eigenvalue problem, Linear Algebra and its Applications, 309 (2000), pp. 339–361.
