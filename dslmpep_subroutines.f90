MODULE dslmpep_subroutines
USE util
IMPLICIT NONE

CONTAINS

!************************************************************************
!			SUBROUTINE DSLM					*
!************************************************************************
! Compute the roots of the scalar polynomial p with real coeffs of	*
! degree d using Laguerre's method. The backward error is stored in be	*
! the number of iterations per root is stored in iter, and the roots are* 
! stored in (er,ei).							*
!************************************************************************
SUBROUTINE dslm(p, er, ei, be, d)
IMPLICIT NONE
!scalar arguments
INTEGER, INTENT(IN) :: d
!array arguments
REAL(dp), INTENT(IN) :: p(*)
REAL(dp), INTENT(INOUT) :: be(*), er(*), ei(*)
!local scalars
LOGICAL :: check
INTEGER :: i, it, nzr, nir, td
REAL(dp) :: tol
!local arrays
REAL(dp), DIMENSION(d+1) :: alpha
!intrinsic procedures
INTRINSIC :: DABS, DBLE, DSQRT, MAXVAL, SUM

!infinite roots and true degree
nir=0
alpha(1:d+1)=DABS(p(1:d+1))
tol=MAXVAL(alpha)
DO i=d+1,2,-1
  IF(alpha(i)<eps*tol) THEN
    er(i-1)=big;ei(i-1)=zero
    nir=nir+1
  ELSE
    EXIT
  ENDIF
ENDDO
td=d-nir
!degree 1 and 2
IF(td==0) THEN
  PRINT*, 'warning: constant polynomial'
ELSEIF(td==1) THEN
  er(1)=-p(1)/p(2)
  ei(1)=zero
ELSEIF(td==2) THEN
  IF(p(2)**2>=4*p(1)*p(3)) THEN
    er(1)=(-p(2)+DSQRT(p(2)**2-4*p(1)*p(3)))/(2*p(3))
    ei(1)=zero
    er(2)=(-p(2)-DSQRT(p(2)**2-4*p(1)*p(3)))/(2*p(3))
    ei(2)=zero
  ELSE
    er(1)=-p(2)/(2*p(3))
    ei(1)=DSQRT(4*p(1)*p(3)-p(2)**2)/(2*p(3))
    er(2)=er(1)
    ei(2)=-ei(1)
  ENDIF
ELSE
  !initial estiamtes
  CALL dsstart(alpha,er,ei,td)
  !zero roots
  nzr=0
  DO i=1,d
    IF(alpha(i)<eps*tol) THEN
      er(i)=zero;ei(i)=zero
      nzr=nzr+1
    ELSE
      EXIT
    ENDIF
  ENDDO
  !backward error for zero and infinite roots
  be(1:nzr)=zero; be(td+1:d)=zero
  !stopping criteria alpha (similar to POLZEROS)
  DO i=1,d+1
    alpha(i)=alpha(i)*(1+3.8*(i-1))
  ENDDO
  !Laguerre's method
  tol=eps*(DSQRT(er(td)**2+ei(td)**2)+1)
  DO i=nzr+1,td
    check=.FALSE.
    DO it=1,itmax
      IF(DABS(er(i))<DSQRT(er(i)**2+ei(i)**2)*eps) THEN
        CALL dslcorr(p, alpha, er, ei, be(i), tol, check, d, i)
      ELSE
        CALL zslcorr(p, alpha, er, ei, be(i), tol, check, d, i)
      ENDIF
      IF(check) THEN
        EXIT
      ENDIF
    ENDDO
  ENDDO
ENDIF
RETURN
END SUBROUTINE dslm

!************************************************************************
!			SUBROUTINE DSLCORR				*
!************************************************************************
! Compute the Laguerre correction term for scalar polynomial p with real*
! coeffs of degree d at real number (tr,0). The magnitude of the	*
! coeffs are stored in alpha, all n previously found roots are stored	*
! in (er,ei). Result is returned in complex number (tr,ti), and backward*
! error of current approximation is stored in be.			*
!************************************************************************
SUBROUTINE dslcorr(p, alpha, er, ei, be, tol, check, d, i)
IMPLICIT NONE
!scalar arguments
LOGICAL, INTENT(INOUT) :: check
INTEGER, INTENT(IN) :: d, i
REAL(dp), INTENT(IN) :: tol
REAL(dp), INTENT(INOUT) :: be
!array arguments
REAL(dp), INTENT(IN) :: p(*), alpha(*)
REAL(dp), INTENT(INOUT) :: er(*), ei(*)
!local scalars
INTEGER :: td
REAL(dp) :: a, b, c, t
COMPLEX(dp) :: x1, x2, y1, y2
!intrinsic procedures
INTRINSIC :: DBLE, DCMPLX, DIMAG, DSQRT, SUM, ZABS, ZSQRT

!initiate variables
t=er(i)
x1=czero; x2=czero
DO td=1,i-1
  y1=DCMPLX(t-er(td),-ei(td))**(-1)
  x1=x1+y1
  x2=x2+y1**2
ENDDO
!split into 2 cases
IF(DABS(t)>1) THEN
  !compute a=revp, be
  CALL drevseval(p, t, a, d, 0)
  CALL drevseval(alpha, DABS(t), be, d, 0)
  be=DABS(a)/be
  IF(be<eps) THEN
    check=.TRUE.
    er(i)=zero
    RETURN
  ELSE
    !solve for b=revp', c=revp''
    CALL drevseval(p, t, b, d, 1)
    CALL drevseval(p, t, c, d, 2)
    !compute y1=p'/p and y2=(p'/p)'
    t=t**(-1)
    y1=t*(d-t*(b/a))
    y2=t**2*(d-2*t*(b/a)+t**2*((b/a)**2-c/a))
  ENDIF
ELSE
  !compute a=p, be
  CALL dseval(p, t, a, d, 0)
  CALL dseval(alpha, DABS(t), be, d, 0)
  be=DABS(a)/be
  IF(be<eps) THEN
    check=.TRUE.
    er(i)=zero
    RETURN
  ELSE
    !solve for b=p', c=p''
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
td=d-(i-1)
!denominator of Laguerre correction term
y1=ZSQRT((td-1)*(td*x2-x1**2))
y2=x1-y1;y1=x1+y1
IF(ZABS(y1)>=ZABS(y2)) THEN
  y1=td*y1**(-1)
  IF(ZABS(y1)<tol) THEN
    check=.TRUE.
    ei(i)=zero
  ELSE
    er(i)=er(i)-DBLE(y1)
    ei(i)=ei(i)-DIMAG(y1)
  ENDIF
ELSE
  y2=td*y2**(-1)
  IF(ZABS(y2)<tol) THEN
    check=.TRUE.
    ei(i)=zero
  ELSE
    er(i)=er(i)-DBLE(y2)
    ei(i)=ei(i)-DIMAG(y2)
  ENDIF
ENDIF
RETURN
END SUBROUTINE dslcorr

!************************************************************************
!			SUBROUTINE ZSLCORR				*
!************************************************************************
! Compute the Laguerre correction term for scalar polynomial p with real*
! coeffs of degree d at complex number (tr,ti). The magnitude of the	*
! coeffs are stored in alpha, all n previously found roots are stored	*
! in (er,ei). Result is return in complex number (tr,ti), and backward	*
! error of current approximation is stored in be.			*
!************************************************************************
SUBROUTINE zslcorr(p, alpha, er, ei, be, tol, check, d, i)
IMPLICIT NONE
!scalar arguments
LOGICAL, INTENT(INOUT) :: check
INTEGER, INTENT(IN) :: d, i
REAL(dp), INTENT(INOUT) :: tol
REAL(dp), INTENT(INOUT) :: be
!array arguments
REAL(dp), INTENT(IN) :: p(*), alpha(*)
REAL(dp), INTENT(INOUT) :: er(*), ei(*)
!local scalars
INTEGER :: td
COMPLEX(dp) :: a, b, c, t, x1, x2, y1, y2
!intrinsic procedures
INTRINSIC :: DBLE, DCMPLX, DIMAG, DSQRT, SUM, ZABS, ZSQRT

!initiate variables
t=DCMPLX(er(i),ei(i))
x1=czero; x2=czero
DO td=1,i-1
  y1=(t-DCMPLX(er(td),ei(td)))**(-1)
  x1=x1+y1
  x2=x2+y1**2
ENDDO
!split into 2 cases
IF(ZABS(t)>1) THEN
  !compute a=revp, be
  CALL zrevseval(p, t, a, d, 0)
  CALL drevseval(alpha, ZABS(t), be, d, 0)
  be=ZABS(a)/be
  IF(be<eps) THEN
    check=.TRUE.
    RETURN
  ELSE
    !solve for b=revp', c=revp''
    CALL zrevseval(p, t, b, d, 1)
    CALL zrevseval(p, t, c, d, 2)
    !compute y1=p'/p and y2=(p'/p)'
    t=t**(-1)
    y1=t*(d-t*(b/a))
    y2=t**2*(d-2*t*(b/a)+t**2*((b/a)**2-c/a))
  ENDIF
ELSE
  !compute a=p, be
  CALL zseval(p, t, a, d, 0)
  CALL dseval(alpha, ZABS(t), be, d, 0)
  be=ZABS(a)/be
  IF(be<eps) THEN
    check=.TRUE.
    RETURN
  ELSE
    !solve for b=p', c=p''
    CALL zseval(p, t, b, d, 1)
    CALL zseval(p, t, c, d, 2)
    !compute y1=p'/p and y2=(p'/p)'
    y1=b/a
    y2=y1**2-(c/a)
  ENDIF
ENDIF
!remove previously found roots
x1=y1-x1
x2=y2-x2
td=d-(i-1)
!denominator of Laguerre correction term
y1=ZSQRT((td-1)*(td*x2-x1**2))
y2=x1-y1;y1=x1+y1
IF(ZABS(y1)>=ZABS(y2)) THEN
  y1=td*y1**(-1)
  IF(ZABS(y1)<tol) THEN
    check=.TRUE.
  ELSE
    er(i)=er(i)-DBLE(y1)
    ei(i)=ei(i)-DIMAG(y1)
  ENDIF
ELSE
  y2=td*y2**(-1)
  IF(ZABS(y2)<tol) THEN
    check=.TRUE.
  ELSE
    er(i)=er(i)-DBLE(y2)
    ei(i)=ei(i)-DIMAG(y2)
  ENDIF
ENDIF
RETURN
END SUBROUTINE zslcorr

!************************************************************************
!			SUBROUTINE DREVSEVAL				*
!************************************************************************
! Evaluate reversal of scalar polynomial p with real coeffs of degree d,*
! and its der=0,1,2 derivatives at real number 1/t, where |t|>1. 	*
! Returns evaluation in a.						*
!************************************************************************
SUBROUTINE drevseval(p, t, a, d, der)
IMPLICIT NONE
!scalar arguments
INTEGER, INTENT(IN) :: d, der
REAL(dp), INTENT(IN) :: t
REAL(dp), INTENT(INOUT) :: a
!array arguments
REAL(dp), INTENT(IN) :: p(*)
!local scalars
INTEGER :: k
REAL(dp) :: y

y=t**(-1)
IF(der==0) THEN
  a=p(1)
  DO k=2,d+1
    a=y*a+p(k)
  ENDDO
ELSEIF(der==1) THEN
  a=d*p(1)
  DO k=2,d
    a=y*a+(d-k+1)*p(k)
  ENDDO
ELSE
  a=d*(d-1)*P(1)
  DO k=2,d-1
    a=y*a+(d-k+1)*(d-k)*p(k)
  ENDDO
ENDIF
RETURN
END SUBROUTINE drevseval


!************************************************************************
!			SUBROUTINE DSEVAL				*
!************************************************************************
! Evaluate scalar polynomial p with real coeffs of degree d, and its	*
! der=0,1,2 derivatives at real number t, where |t|<=1. Returns		*
! evaluation in a.							*
!************************************************************************
SUBROUTINE dseval(p, t, a, d, der)
IMPLICIT NONE
!scalar arguments
INTEGER, INTENT(IN) :: d, der
REAL(dp), INTENT(IN) :: t
REAL(dp), INTENT(INOUT) :: a
!array arguments
REAL(dp), INTENT(IN) :: p(*)
!local scalars
INTEGER :: k

IF(der==0) THEN
  a=p(d+1)
  DO k=d,1,-1
    a=t*a+p(k)
  ENDDO
ELSEIF(der==1) THEN
  a=d*p(d+1)
  DO k=d,2,-1
    a=t*a+(k-1)*p(k)
  ENDDO
ELSE
  a=d*(d-1)*p(d+1)
  DO k=d,3,-1
    a=t*a+(k-1)*(k-2)*p(k)
  ENDDO
ENDIF
RETURN
END SUBROUTINE dseval

!************************************************************************
!			SUBROUTINE ZREVSEVAL				*
!************************************************************************
! Evaluate reversal of scalar polynomial p with real coeffs of degree d,*
! and its der=0,1,2 derivatives at complex number 1/t, where |t|>1. 	*
! Returns evaluation in a.						*
!************************************************************************
SUBROUTINE zrevseval(p, t, a, d, der)
IMPLICIT NONE
!scalar arguments
INTEGER, INTENT(IN) :: d, der
COMPLEX(dp), INTENT(IN) :: t
COMPLEX(dp), INTENT(INOUT) :: a
!array arguments
REAL(dp), INTENT(IN) :: p(*)
!local scalars
INTEGER :: k
COMPLEX(dp) :: y

y=t**(-1)
IF(der==0) THEN
  a=p(1)
  DO k=2,d+1
    a=y*a+p(k)
  ENDDO
ELSEIF(der==1) THEN
  a=d*p(1)
  DO k=2,d
    a=y*a+(d-k+1)*p(k)
  ENDDO
ELSE
  a=d*(d-1)*P(1)
  DO k=2,d-1
    a=y*a+(d-k+1)*(d-k)*p(k)
  ENDDO
ENDIF
RETURN
END SUBROUTINE zrevseval

!************************************************************************
!			SUBROUTINE ZSEVAL				*
!************************************************************************
! Evaluate scalar polynomial p with real coeffs of degree d, and its	*
! der=0,1,2 derivatives at complex number t, where |t|<=1. Returns 	*
! evaluation in a.							*
!************************************************************************
SUBROUTINE zseval(p, t, a, d, der)
IMPLICIT NONE
!scalar arguments
INTEGER, INTENT(IN) :: d, der
COMPLEX(dp), INTENT(IN) :: t
COMPLEX(dp), INTENT(INOUT) :: a
!array arguments
REAL(dp), INTENT(IN) :: p(*)
!local scalars
INTEGER :: k

IF(der==0) THEN
  a=p(d+1)
  DO k=d,1,-1
    a=t*a+p(k)
  ENDDO
ELSEIF(der==1) THEN
  a=d*p(d+1)
  DO k=d,2,-1
    a=t*a+(k-1)*p(k)
  ENDDO
ELSE
  a=d*(d-1)*p(d+1)
  DO k=d,3,-1
    a=t*a+(k-1)*(k-2)*p(k)
  ENDDO
ENDIF
RETURN
END SUBROUTINE zseval

!************************************************************************
!			SUBROUTINE DSSTART				*
!************************************************************************
! Compute the initial esimates of the roots of a real scalar polynomial *
! using the Newton Polygon method. The absolute value of the 		*
! coefficients of the polynomial of degree d are stored in alpha, the 	*
! real and imaginary parts of the initial estiamtes are returned in er 	*
! and ei.         							*
!************************************************************************
SUBROUTINE dsstart(alpha, er, ei, d)
IMPLICIT NONE
!scalar arguments
INTEGER, INTENT(IN) :: d
!array arguments
REAL(dp), INTENT(IN) :: alpha(*)
REAL(dp), INTENT(INOUT) :: er(*), ei(*)
!parameters
REAL(dp), PARAMETER :: pi2 = 6.2831853071795865_dp, sigma = 0.7_dp
!local scalars
INTEGER :: c, i, iold, j, nzeros
REAL(dp) :: ang, r, th
!local arrays
LOGICAL, DIMENSION(d+1) :: h
REAL(dp), DIMENSION(d+1) :: a
!intrinsic procedures
INTRINSIC :: DCOS, DEXP, DLOG, DSIN

!compute log(alpha)
DO I=1,d+1
  IF(alpha(i)>=eps) THEN
    a(i)=DLOG(alpha(i))
  ELSE
    a(i)=-one
  ENDIF
ENDDO
!compute upper convex hull
CALL cnvex(d+1,a,h)
!compute initial estimates
iold=1;c=0;th=pi2/d
DO i=2,d+1
  IF(h(i)) THEN
    nzeros=i-iold
    r=DEXP((a(iold)-a(i))/nzeros)
    ang=pi2/nzeros
    DO j=1,nzeros
      er(c+j)=r*DCOS(ang*j+th*i+sigma)
      ei(c+j)=r*DSIN(ang*j+th*i+sigma)
    ENDDO
    c=c+nzeros
    iold=i
  ENDIF
ENDDO
RETURN
END SUBROUTINE dsstart

END MODULE dslmpep_subroutines
