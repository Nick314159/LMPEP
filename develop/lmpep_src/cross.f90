!>\author Thomas R. Cameron, Davidson College
!>\author Nikolas I. Steckley, Portland State University
!>\date 2017
!>\brief <b> CROSS computes the angle between three points. </b>
!>\par Purpose:
!>\verbatim
!> CROSS uses the determinant of a 2x2 matrix.
!>\endverbatim
!>\param[in] a
!>\verbatim Double precision array of dimension n, contains the log of the moduli of polynomial coefficients,
!> ordered from constant to leading. \endverbatim
!>\param[in] n
!>\verbatim  Integer, size of a.\endverbatim
!>\param[in] h
!>\verbatim  Integer array that contains the indexes corresponding to the upper envelope of the convex hull.\endverbatim
!>\param[in] c
!>\verbatim Integer, number of indeces in the upper envelop of the convex hull.\endverbatim
!>\param[in] i
!>\verbatim Integer, current index to check.\endverbatim
!************************************************************************
DOUBLE PRECISION FUNCTION cross(h, a, c, i)
IMPLICIT NONE
! scalar arguments
INTEGER, INTENT(IN) 			:: c, i
! array arguments
INTEGER, INTENT(IN) 			:: h(*)
DOUBLE PRECISION, INTENT(IN)	:: a(*)

cross = (a(i)-a(h(c-1)))*(h(c)-h(c-1)) - (a(h(c))-a(h(c-1)))*(i-h(c-1))
RETURN
END FUNCTION cross