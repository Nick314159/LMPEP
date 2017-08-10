SUBROUTINE dslcorr(p, alpha, tol, check, deg, ind, er, ei, berr)
USE util
IMPLICIT NONE
!scalar arguments
LOGICAL, INTENT(INOUT) :: check
INTEGER(KIND=in1), INTENT(IN) :: deg
INTEGER, INTENT(IN) :: d, i
REAL(KIND=re8), INTENT(IN) :: tol
REAL(KIND=re8), INTENT(INOUT) :: berr
!array arguments
REAL(KIND=re8), INTENT(IN) :: p(*), alpha(*)
REAL(KIND=re8), INTENT(INOUT) :: er(*), ei(*)
!local scalars
INTEGER :: k
REAL(dp) :: a, b, c, t
DOUBLE COMPLEX :: x1, x2, y1, y2
!intrinsic procedures
INTRINSIC :: DABS, DBLE, DCMPLX, DIMAG, ZABS, ZSQRT

!initiate variables
x1=czero; x2=czero
DO k=1,i-1
  y1=1/DCMPLX(er(i)-er(k),-ei(k))
  x1=x1+y1
  x2=x2+y1**2
ENDDO
t=er(i)
!split into 2 cases
IF(DABS(t)>1) THEN
  !compute a=revp, berr
  t=1/t
  CALL drevseval(p, t, a, d, 0)
  CALL drevseval(alpha, DABS(t), berr, d, 0)
  berr=MIN(DABS(a)/berr, DABS(a))
  IF(berr<eps) THEN
    ei(i)=zero
    check=.TRUE.
    RETURN
  ELSE
    !compute b=revp', c=revp''
    CALL drevseval(p, t, b, d, 1)
    CALL drevseval(p, t, c, d, 2)
    !compute y1=p'/p and y2=(p'/p)'
    y1=t*(d-t*(b/a))
    y2=t**2*(d-2*t*(b/a)+t**2*((b/a)**2-c/a))
  ENDIF
ELSE
  !compute a=p, berr
  CALL dseval(p, t, a, d, 0)
  CALL dseval(alpha, DABS(t), berr, d, 0)
  berr=MIN(DABS(a)/berr, DABS(a))
  IF(berr<eps) THEN
    ei(i)=zero
    check=.TRUE.
    RETURN
  ELSE
    !compute b=p', c=p''
    CALL dseval(p, t, b, d, 1)
    CALL dseval(p, t, c, d, 2)
    !compute y1=p'/p and y2=(p'/p)'
    y1=b/a
    y2=y1**2-c/a
  ENDIF
ENDIF
!remove previously found roots
x1=y1-x1
x2=y2-x2
k=d-i+1
!denominator of Laguerre correction term
y1=ZSQRT((k-1)*(k*x2-x1**2))
y2=x1-y1; y1=x1+y1
IF(ZABS(y1)>=ZABS(y2)) THEN
  y1=k/y1
  IF(ZABS(y1)<tol) THEN
    ei(i)=zero
    check=.TRUE.
  ELSE
    er(i)=er(i)-DBLE(y1)
    ei(i)=-DIMAG(y1)
  ENDIF
ELSE
  y2=k/y2
  IF(ZABS(y2)<tol) THEN
    ei(i)=zero
    check=.TRUE.
  ELSE
    er(i)=er(i)-DBLE(y2)
    ei(i)=-DIMAG(y2)
  ENDIF
ENDIF
RETURN
END SUBROUTINE dslcorr
