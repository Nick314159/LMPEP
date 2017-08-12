!>\author Thomas R. Cameron, Davidson College
!>\author Nikolas I. Steckley, Portland State University
!>\date 2017
!>\brief <b> DSLCORR computes the Lagueree correction term of real polynomial with real root approximation, while ``dividing out'' previously computed roots. </b>
!>\par Purpose:
!>\verbatim
!> DSLCORR calculates y1=p'(t)/p(t) and y2=-(p'(t)/p(t))', where t=er(ind) is the current eigenvalue approximation, and then updates y1 and y2 by subtracting the appropriate sum of preivously computed roots indexed by 1,...ind-1. From there the Laguerre correction term can be computed as deg/(y1+-sqrt((deg-1)*(deg*y2-y1**2)), where +- is choosen to maximize the denominator. 
!>\endverbatim
!>\param[in] p
!>\verbatim Double precision array of dimension (deg+1), contains polynomial coefficients, ordered from constant to leading. \endverbatim
!>\param[in] alpha
!>\verbatim Double precision array of dimension (deg+1), contains moduli of polynomial coefficients.\endverbatim
!>\param[in] tol
!>\verbatim Double precision, used to determine potential convergence of eigenvalue approximation..\endverbatim
!>\param[in,out] conv
!>\verbatim Logical, returns true if convergence is reached.\endverbatim
!>\param[in] deg
!>\verbatim  Integer, degree of the polynomial.\endverbatim
!>\param[in] ind
!>\verbatim  Integer, index of current eigenvalue approximation.\endverbatim
!>\param[in,out] er
!>\verbatim  Double precision array of dimension deg, real part of eigenvalue approximations.\endverbatim
!>\param[in,out] ei
!>\verbatim  Double precision array of dimension deg, imaginary part of eigenvalue approximations.\endverbatim
!>\param[in,out] berr
!>\verbatim  Double precision number, backward error in current eigenvalue approximation.\endverbatim
!>\note MEMORY: O(deg), FLOPS: O(deg)
!************************************************************************
SUBROUTINE dslcorr(p, alpha, tol, conv, deg, ind, er, ei, berr)
IMPLICIT NONE
!scalar arguments
LOGICAL, INTENT(INOUT)          :: conv
INTEGER, INTENT(IN)             :: deg, ind
DOUBLE PRECISION, INTENT(IN)    :: tol
DOUBLE PRECISION, INTENT(INOUT) :: berr
!array arguments
DOUBLE PRECISION, INTENT(IN)    :: p(*), alpha(*)
DOUBLE PRECISION, INTENT(INOUT) :: er(*), ei(*)
!parameters
DOUBLE PRECISION, PARAMETER     :: zero=0.0D0
DOUBLE PRECISION, PARAMETER     :: eps=EPSILON(zero)
!local scalars
INTEGER                         :: k
DOUBLE PRECISION                :: a, b, c, t
DOUBLE COMPLEX                  :: x1, x2, y1, y2
!intrinsic procedures
INTRINSIC                       :: DABS, DBLE, DCMPLX, DIMAG, ZABS, ZSQRT
!external subroutines
EXTERNAL                        :: dseval, drevseval

!initiate variables
x1=zero; x2=zero
DO k=1,ind-1
  y1=DCMPLX(er(ind)-er(k),-ei(k))**(-1)
  x1=x1+y1
  x2=x2+y1**2
ENDDO
t=er(ind)
!split into 2 cases
IF(DABS(t)>1) THEN
  !compute a=revp, berr
  t=t**(-1)
  CALL drevseval(p, t, deg, 0, a)
  CALL drevseval(alpha,DABS(t),deg,0,berr)
  berr=DABS(a)*berr**(-1)
  IF(berr<eps) THEN
    ei(ind)=zero
    conv=.TRUE.
    RETURN
  ENDIF
  !compute b=revp', c=revp''
  CALL drevseval(p, t, deg, 1, b)
  CALL drevseval(p, t, deg, 2, c)
  !compute y1=p'/p and y2=(p'/p)'
  y1=t*(deg-t*(b*a**(-1)))
  y2=t**2*(deg-2*t*(b*a**(-1))+t**2*((b*a**(-1))**2-(c*a**(-1))))
ELSE
  !compute a=p, berr
  CALL dseval(p, t, deg, 0, a)
  CALL dseval(alpha, DABS(t), deg, 0, berr)
  berr=MIN(DABS(a)*berr**(-1), DABS(a))
  IF(berr<eps) THEN
    ei(ind)=zero
    conv=.TRUE.
    RETURN
  ENDIF
  !compute b=p', c=p''
  CALL dseval(p, t, deg, 1, b)
  CALL dseval(p, t, deg, 2, c)
  !compute y1=p'/p and y2=(p'/p)'
  y1=b*a**(-1)
  y2=y1**2-(c*a**(-1))
ENDIF
!remove previously found roots
x1=y1-x1
x2=y2-x2
!denominator of Laguerre correction term
y1=ZSQRT((deg-1)*(deg*x2-x1**2))
y2=x1-y1
y1=x1+y1
!choose term that maximizes denominator
IF(ZABS(y1)>ZABS(y2)) THEN
  y1=deg*y1**(-1)
  IF(ZABS(y1)<tol) THEN
    ei(ind)=zero
    conv=.TRUE.
  ELSE
    er(ind)=er(ind)-DBLE(y1)
    ei(ind)=-DIMAG(y1)
  ENDIF
ELSE
  y2=deg*y2**(-1)
  IF(ZABS(y2)<tol) THEN
    ei(ind)=zero
    conv=.TRUE.
  ELSE
    er(ind)=er(ind)-DBLE(y2)
    ei(ind)=-DIMAG(y2)
  ENDIF
ENDIF
RETURN
END SUBROUTINE dslcorr
