!>\author Thomas R. Cameron, Davidson College
!>\author Nikolas I. Steckley, Portland State University
!>\date 2017
!>\brief <b> DSSTART computes initial estimates to the roots of a polynomial. </b>
!>\par Purpose:
!>\verbatim
!> DSSTART uses the Newton Polygon of a polynomial to compute initial estimates to that polynomials roots.  
!>\endverbatim
!>\param[in] alpha
!>\verbatim Double precision array of dimension (deg+1), contains moduli of polynomial coefficients,
!> ordered from constant to leading. \endverbatim
!>\param[in] deg
!>\verbatim  Integer, degree of the polynomial.\endverbatim
!>\param[out] er
!>\verbatim  Double precision array of dimension deg, real part of eigenvalue approximations.\endverbatim
!>\param[out] ei
!>\verbatim  Double precision array of dimension deg, imaginary part of eigenvalue approximations.\endverbatim
!************************************************************************
SUBROUTINE dsstart(alpha, deg, er, ei)
IMPLICIT NONE
! scalar arguments
INTEGER, INTENT(IN)                 :: deg
! array arguments
DOUBLE PRECISION, INTENT(IN)        :: alpha(*)
DOUBLE PRECISION, INTENT(OUT)       :: er(*), ei(*)
! local scalars
INTEGER                             :: c, i, j, k, nzeros
DOUBLE PRECISION                    :: ang, r, th
! local arrays
INTEGER, DIMENSION(deg+1)           :: h
DOUBLE PRECISION, DIMENSION(deg+1)  :: a
! intrinsic functions
INTRINSIC                           :: dlog, dcos, dsin
! external subroutines
EXTERNAL							:: conv_hull 
! parameters
DOUBLE PRECISION, PARAMETER         :: zero = 0.0D0, one = 1.0D0, big_one = -1.0D30
DOUBLE PRECISION, PARAMETER         :: pi2 = 6.283185307179586D0, sigma = 0.7D0

! compute log(alpha)
DO i=1,deg+1
	IF(alpha(i)>zero) THEN
		a(i) = dlog(alpha(i))
	ELSE
		a(i) = big_one
	ENDIF
ENDDO
! compute upper convex hull
CALL conv_hull(deg+1, a, h, c)
! compute initial estimates
k=0; th=pi2/deg
DO i=c-1,1,-1
	nzeros = h(i)-h(i+1)
	r = (alpha(h(i+1))/alpha(h(i)))**(one/nzeros)
	ang = pi2/nzeros
	DO j=1,nzeros
		er(k+j) = r*dcos(ang*j+th*h(i)+sigma)
		ei(k+j) = r*dsin(ang*j+th*h(i)+sigma)
	ENDDO
	k = k+nzeros
ENDDO
RETURN
END SUBROUTINE dsstart