SUBROUTINE dsacorr(p, alpha, tol, deg, ind, check, er, ei, berr)
IMPLICIT NONE
!scalar arguments
LOGICAL, INTENT(OUT)            :: check
INTEGER, INTENT(IN)             :: deg, ind
DOUBLE PRECISION, INTENT(IN)    :: tol
DOUBLE PRECISION, INTENT(OUT)   :: berr
!array arguments
DOUBLE PRECISION, INTENT(IN)    :: p(*), alpha(*)
DOUBLE PRECISION, INTENT(INOUT) :: er(*), ei(*)
!local scalars
INTEGER                         :: k
DOUBLE PRECISION                :: a, b, t, t2
DOUBLE COMPLEX                  :: x1, y1
!intrinsic procedures
INTRINSIC                       :: dabs, dble, dcmplx, dimag, epsilon, zabs, zsqrt
!parameters
DOUBLE PRECISION, PARAMETER     :: zero=0.0D0
DOUBLE PRECISION, PARAMETER     :: eps=epsilon(zero)
DOUBLE COMPLEX, PARAMETER       :: czero=dcmplx(zero,zero)
!external subroutines
EXTERNAL                        :: dseval, drevseval

!initiate variables
x1=czero
DO k=1,ind-1
  y1=dcmplx(er(ind)-er(k),-ei(k))**(-1)
  x1=x1+y1
ENDDO
DO k=ind+1,deg
  y1=dcmplx(er(ind)-er(k),-ei(k))**(-1)
  x1=x1+y1
ENDDO
t=er(ind)
t2=dabs(t)
!split into 2 cases
IF(t2>1) THEN
  !compute a=revp, berr
  t=t**(-1)
  t2=t2**(-1)
  CALL drevseval(p, t, deg, 0, a)
  CALL drevseval(alpha, t2, deg, 0, berr)
  berr=dabs(a)*berr**(-1)
  IF(berr<eps) THEN
    ei(ind)=zero
    check=.FALSE.
    RETURN
  ENDIF
  !compute b=revp'
  CALL drevseval(p, t, deg, 1, b)
  !compute y1=p'/p
  b=b*a**(-1)
  y1=t*(deg-t*b)
ELSE
  !compute a=p, berr
  CALL dseval(p, t, deg, 0, a)
  CALL dseval(alpha, t2, deg, 0, berr)
  berr=dabs(a)*berr**(-1)
  IF(berr<eps) THEN
    ei(ind)=zero
    check=.FALSE.
    RETURN
  ENDIF
  !compute b=p'
  CALL dseval(p, t, deg, 1, b)
  !compute y1=p'/p
  y1=b*a**(-1)
ENDIF
!remove previously found roots (Aberth correction term)
x1=y1-x1
y1=x1**(-1)
!Aberth iterate
er(ind)=er(ind)-dble(y1)
ei(ind)=-dimag(y1)
RETURN
END SUBROUTINE dsacorr
