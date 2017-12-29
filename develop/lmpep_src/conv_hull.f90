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
!************************************************************************
SUBROUTINE conv_hull(n, a, h, c)
IMPLICIT NONE
! scalar argument
INTEGER, INTENT(IN) 			:: n
INTEGER, INTENT(INOUT) 			:: c
! array argument
INTEGER, INTENT(INOUT)			:: h(*)
DOUBLE PRECISION, INTENT(IN)	:: a(*)
! local scalars
INTEGER 						:: i
! parameters
DOUBLE PRECISION, PARAMETER 	:: toler=0.0D0
! external functions
DOUBLE PRECISION 				:: cross
EXTERNAL 						:: cross

! build upper envelop of convex hull
 c=0
DO i=n,1,-1
	DO WHILE(c>=2 .and. cross(h,a,c,i)<toler)
		c = c-1
	END DO
	c = c+1
	h(c) = i
END DO
RETURN
END SUBROUTINE conv_hull