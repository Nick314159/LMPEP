SUBROUTINE dzslcorr2(p, alpha, tol, deg, ind, check, er, ei, berr)
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
DOUBLE PRECISION                :: t2
DOUBLE COMPLEX                  :: a, b, c, t, x1, x2, y1, y2
!intrinsic procedures
INTRINSIC                       :: dabs, dble, dcmplx, dimag, epsilon, zabs, zsqrt
!parameters
DOUBLE PRECISION, PARAMETER     :: eps=epsilon(0.0D0)
DOUBLE COMPLEX, PARAMETER       :: czero=dcmplx(0.0D0,0.0D0)
!external subroutines
EXTERNAL                        :: dzseval, dzrevseval
!external functions
DOUBLE PRECISION                :: dzmod
EXTERNAL                        :: dzmod

!initiate variables
x1=czero; x2=czero
DO k=1,ind-1
  y1=dcmplx(er(ind)-er(k),ei(ind)-ei(k))**(-1)
  x1=x1+y1
  x2=x2+y1**2
ENDDO
DO k=ind+1,deg
  y1=dcmplx(er(ind)-er(k),ei(ind)-ei(k))**(-1)
  x1=x1+y1
  x2=x2+y1**2
ENDDO
x2=deg*x2-deg*x1**2*(deg-1)**(-1)
t=dcmplx(er(ind),ei(ind))
t2=dzmod(er(ind),ei(ind))
!split into 2 cases
IF(t2>1) THEN
  !compute a=revp, berr
  t=t**(-1)
  t2=t2**(-1)
  CALL dzrevseval(p, t, deg, 0, a)
  CALL drevseval(alpha, t2, deg, 0, berr)
  berr=zabs(a)*berr**(-1)
  IF(berr<eps) THEN
    check=.FALSE.
    RETURN
  ENDIF
  !compute b=revp', c=revp''
  CALL dzrevseval(p, t, deg, 1, b)
  CALL dzrevseval(p, t, deg, 2, c)
  !compute y1=p'/p and y2=-(p'/p)'
  b=b*a**(-1)
  c=c*a**(-1)
  y1=t*(deg-t*b)
  y2=t**2*(deg-2*t*b+t**2*(b**2-c))
ELSE
  !compute a=p, berr
  CALL dzseval(p, t, deg, 0, a)
  CALL dseval(alpha, t2, deg, 0, berr)
  berr=zabs(a)*berr**(-1)
  IF(berr<eps) THEN
    check=.FALSE.
    RETURN
  ENDIF
  !compute b=p', c=p''
  CALL dzseval(p, t, deg, 1, b)
  CALL dzseval(p, t, deg, 2, c)
  !compute y1=p'/p and y2=-(p'/p)'
  b=b*a**(-1)
  c=c*a**(-1)
  y1=b
  y2=b**2-c
ENDIF
!denominator of Laguerre correction term
x1=zsqrt((deg-1)*(deg*y2-y1**2-x2))
y2=y1-x1
y1=y1+x1
!choose term that maximizes denominator
IF(zabs(y1)>zabs(y2)) THEN
  y1=deg*y1**(-1)
  IF(zabs(y1)<tol) THEN
    check=.FALSE.
  ELSE
    er(ind)=er(ind)-dble(y1)
    ei(ind)=ei(ind)-dimag(y1)
  ENDIF
ELSE
  y2=deg*y2**(-1)
  IF(zabs(y2)<tol) THEN
    check=.FALSE.
  ELSE
    er(ind)=er(ind)-dble(y2)
    ei(ind)=ei(ind)-dimag(y2)
  ENDIF
ENDIF
RETURN
END SUBROUTINE dzslcorr2
