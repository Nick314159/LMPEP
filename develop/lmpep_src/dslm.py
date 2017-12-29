"""
author Thomas R. Cameron, Davidson College
author Nikolas I. Steckley, Portland State University
date 2017
DSLM simultaneously computes the roots of a real polynomial. </b>
Purpose: DSLM computes the roots of a real polynomial simultaneously using a modified Laguerre's method. 
endverbatim
param[in] p
verbatim Double precision array of dimension (deg+1), contains polynomial coefficients, ordered from constant to leading. \endverbatim
param[in] deg
verbatim Integer, degree of the polynomial.\endverbatim
param[out] er
verbatim Double precision array of dimension deg, real part of eigenvalue approximations.\endverbatim
param[out] ei
verbatim  Double precision array of dimension deg, imaginary part of eigenvalue approximations.\endverbatim
param[out] berr
verbatim  Double precision array of dimension deg, backward error of each eigenvalue approximation.\endverbatim 
"""
import sys
itmax = 50
epsilon = sys.float_info.epsilon

def dslm(p):
#def dslm(p, deg, er, ei, berr):
    deg = len(p)
    #Alpha
    alpha = [abs(p[i]) for i in xrange(1, deg+1)]
    #Checks
    checks = [True for i in xrange(1, deg)]
    #Initial Estimates
    er, ei = dsstart(alpha)
    #Lagurre's method
    for it in xrange(1, itmax):
        for i in xrange(1, deg):
            if check[i]:
                tol = epislon * dzmod(er([i], ei[i])
                if abs(ei[i] < tol:
                   dslcorr(p, alpha, tol, deg, i, check(i), er, ei, berr(i))
                else:
                    dzslcorr(p, alpha, tol, deg, i, check(i), er, ei, berr(i))
                
            
        
    


