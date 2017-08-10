!************************************************************************
!                           SUBROUTINE DSEVAL				            *
!           Authors: Thomas R. Cameron, Nikolas I. Steckley             *
!                           Date: 8/10/2017                             *
!************************************************************************
! Compute Laguree Correction term of real scalar polynomial at real     *
! number while ''dividing out'' previously computed roots.              *
!************************************************************************
! Input Variables:                                                      *
!   p: REAL(re8), array of dimension (deg+1),                           *
!       contains polynomial coefficients from constant to leading.      *
!   alpha: REAL(re8), array of dimension (d+1),                         *
!       contains modulie of polynomial coefficients.                    *
!   tol: REAL(re8), convergence tolerance.                              *
!   conv: LOGICAL, convergence parameter.                               *
!   deg: INTEGER(in4), degree of the polynomial                         *
!   ind: INTEGER(in4), index of eigenvalue being updated                *
!                                                                       *
! Output Variables:                                                     *
!   er: REAL(re8), array of dimension deg,                              *
!       contains the real part of the eigenvalue approximations.        *
!   ei: REAL(re8), array of dimension deg,                              *
!       contains the imaginary part of the eigenvalue approximations.   *
!   berr: REAL(re8), backward error in current eigenvalue approximation.*
!                                                                       *
!   MEMORY: O(deg), FLOPS: O(deg)                                       *
!************************************************************************
SUBROUTINE dslcorr(p, alpha, tol, conv, deg, ind, er, ei, berr)
USE util
IMPLICIT NONE
!scalar arguments
LOGICAL, INTENT(INOUT)          :: conv
INTEGER(KIND=in4), INTENT(IN)   :: deg, ind
REAL(KIND=re8), INTENT(IN)      :: tol
REAL(KIND=re8), INTENT(INOUT)   :: berr
!array arguments
REAL(KIND=re8), INTENT(IN)      :: p(*), alpha(*)
REAL(KIND=re8), INTENT(INOUT)   :: er(*), ei(*)
!local scalars
INTEGER(KIND=in4)               :: k
REAL(KIND=re8)                  :: a, b, c, t
COMPLEX(KIND=re8)               :: x1, x2, y1, y2
!intrinsic procedures
INTRINSIC                       :: DABS, DBLE, DCMPLX, DIMAG, MIN, ZABS, ZSQRT
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
  berr=MIN(DABS(a)*berr**(-1), DABS(a))
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
k=deg-(ind-1)
!denominator of Laguerre correction term
y1=ZSQRT((k-1)*(k*x2-x1**2))
y2=x1-y1
y1=x1+y1
IF(ZABS(y1)>ZABS(y2)) THEN
  y1=k*y1**(-1)
  IF(ZABS(y1)<tol) THEN
    ei(ind)=zero
    conv=.TRUE.
  ELSE
    er(ind)=er(ind)-DBLE(y1)
    ei(ind)=-DIMAG(y1)
  ENDIF
ELSE
  y2=k*y2**(-1)
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
