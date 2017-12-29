"""
author Thomas R. Cameron, Davidson College
author Nikolas I. Steckley, Portland State University
date 2017
brief <b> DSSTART2 computes initial estimates to the roots of a polynomial. </b>
par Purpose:
verbatim
!> DSSTART2 uses the Newton Polygon of a polynomial to compute initial estimates to that polynomials roots.  
endverbatim
param[in] alpha
verbatim Double precision array of dimension (deg+1), contains moduli of polynomial coefficients,
!> ordered from constant to leading. \endverbatim
param[in] deg
verbatim  Integer, degree of the polynomial.\endverbatim
param[out] er
verbatim  Double precision array of dimension deg, real part of eigenvalue approximations.\endverbatim
param[out] ei
verbatim  Double precision array of dimension deg, imaginary part of eigenvalue approximations.\endverbatim 
"""
import math

big_one = -1.0e30
pi2 = pi * 2
one = 1.0
sigma = 0.7
def dsstart(alpha):
    deg = len(alpha) - 1
    a = []
    #Compute log(alpha)
    for i in xrange(1, deg+1):
        if alpha[i] > 0.0:
            a[i] = log(alpha[i])
        else:
            a[i] = big_one
    #compute upper convex hull
    h, c = conv_hull(a)
    #compute initial estimates
    k = 0
    th = pi2 / deg
    for i in xrange(c-1, 1, -1):
        nzeros = h[i] - h[i+1]
        r = alpha[h[i+1]] / alpha[h[i]]**(one/nzeros)
        ang = pi2/nzeros
        for j in xrange(1, nzeros);
            er[k+j] = r * cos(ang*j + th*h[i] + sigma) 
            ei[k+j] = r * sin(ang*j + th*h[i] + sigma)
        k = k + nzeros
    return er, ei


