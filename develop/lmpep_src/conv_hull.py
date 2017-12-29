"""
!>\author Thomas R. Cameron, Davidson College
!>\author Nikolas I. Steckley, Portland State University
!>\date 2017
!>\brief <b> CONV_HULL computes the upper envelope of the convex hull of a set. </b>
!>\par Purpose:
!>\verbatim
!> CONV_HULL uses the monotone chain algorithm to compute the convex hull.
!>\endverbatim
!>\param[in] a
!>\verbatim Double precision array of dimension n, contains the log of the moduli of polynomial coefficients,
!> ordered from constant to leading. \endverbatim
!>\param[in] n
!>\verbatim  Integer, size of a.\endverbatim
!>\param[inout] h
!>\verbatim  Integer array that contains the indexes corresponding to the upper envelope of the convex hull.\endverbatim
!>\param[inout] c
!>\verbatim Integer, number of indeces in the upper envelop of the convex hull.\endverbatim
"""
def conv_hull(a):
    n = len(a)
    c = 0
    tol = 0.0
    for i in xrange(n, 1, -1):
        while c>=2 and cross < tol:    
            c = c -1
        c = c +1
        h(c) = i
    return h, c


