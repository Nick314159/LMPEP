!>\author Thomas R. Cameron, Davidson College
!>\author Nikolas I. Steckley, Portland State University
!>\date 2017
!>\brief <b> DSSTART computes initial estimates to the roots of a polynomial. </b>
!>\par Purpose:
!>\verbatim
!> DSSTART uses the Newton Polygon of a polynomial to compute initial estimates to that polynomials roots.  
!>\endverbatim
!>\param[in] alpha
!>\verbatim Double precision array of dimension (deg+1), contains moduli of polynomial coefficients, ordered from constant to leading. \endverbatim
!>\param[in] deg
!>\verbatim  Integer, degree of the polynomial.\endverbatim
!>\param[in,out] er
!>\verbatim  Double precision array of dimension deg, real part of eigenvalue approximations.\endverbatim
!>\param[out,out] ei
!>\verbatim  Double precision array of dimension deg, imaginary part of eigenvalue approximations.\endverbatim
!>\note MEMORY: O(deg), FLOPS: O(deg)
!************************************************************************
SUBROUTINE dsstart(alpha, deg, er, ei)
IMPLICIT NONE
!scalar arguments
INTEGER, INTENT(IN)                 :: deg
!array arguments
DOUBLE PRECISION, INTENT(IN)        :: alpha(*)
DOUBLE PRECISION, INTENT(INOUT)     :: er(*), ei(*)
!local scalars
INTEGER                             :: c, i, iold, j, nzeros
DOUBLE PRECISION                    :: ang, r, th
!local arrays
LOGICAL, DIMENSION(deg+1)           :: h
DOUBLE PRECISION, DIMENSION(deg+1)  :: a
!intrinsic procedures
INTRINSIC                           :: dcos, dexp, dlog, dsin, epsilon
!parameters
DOUBLE PRECISION, PARAMETER         :: zero=0.0D0, one=1.0D0
DOUBLE PRECISION, PARAMETER         :: eps=epsilon(zero)
DOUBLE PRECISION, PARAMETER         :: pi2 = 6.283185307179586D0, sigma = 0.7D0
!external subroutines
EXTERNAL                            :: cnvex

!compute log(alpha)
DO i=1,deg+1
  IF(alpha(i)>=eps) THEN
    a(i)=dlog(alpha(i))
  ELSE
    a(i)=-one
  ENDIF
ENDDO
!compute upper convex hull
CALL cnvex(deg+1,a,h)
!compute initial estimates
iold=1; c=0; th=pi2/deg
DO i=2,deg+1
  IF(h(i)) THEN
    nzeros=i-iold
    r=dexp((a(iold)-a(i))/nzeros)
    ang=pi2/nzeros
    DO j=1,nzeros
      er(c+j)=r*dcos(ang*j+th*i+sigma)
      ei(c+j)=r*dsin(ang*j+th*i+sigma)
    ENDDO
    c=c+nzeros
    iold=i
  ENDIF
ENDDO
RETURN
END SUBROUTINE dsstart
