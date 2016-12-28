MODULE dgeeam_subroutines
USE dgelmpep_subroutines
IMPLICIT NONE

CONTAINS

!************************************************************************
!			SUBROUTINE DGEEAM				*
!************************************************************************
! Compute the eigenvalues and eigenvectors of the matrix polynomial p 	*
! with real coefficients of degree d and size n using Ehrlich-Aberth	*
! method. The norm of the coefficients of p are stored in ncoeff, iseed *
! is used when computing an upper bound on the backward error, stored in* 
! berr of each eigenpair. Eigenvalues are stored in (er,ei) and 	*
! eigenvectors in (xr,xi) and (yr,yi).					*
!************************************************************************
SUBROUTINE dgeeam(p, xr, xi, yr, yi, er, ei, berr, ncoeff, iseed, d, n, opt)
IMPLICIT NONE
!scalar arguments
CHARACTER(LEN=2), INTENT(IN) :: opt
INTEGER, INTENT(IN) :: d, n
!array arguments
INTEGER, INTENT(INOUT) :: iseed(*)
REAL(dp), INTENT(IN) :: ncoeff(*), p(n,*)
REAL(dp), INTENT(INOUT) :: berr(*), er(*), ei(*)
REAL(dp), INTENT(INOUT) :: xr(n,*), xi(n,*), yr(n,*), yi(n,*)
!local scalars
INTEGER :: c, die, dze, i, it, lwork
REAL(dp) :: tol
!local arrays
LOGICAL, DIMENSION(n*d) :: check
!intrinsic procedures
INTRINSIC :: DABS, MAX

!initalize variables
check=.TRUE.
!initial estimates
CALL dgestart(p, xr, xi, yr, yi, er, ei, ncoeff, d, die, dze, lwork, n, opt)
 c=die+dze
!Ehrlich-Aberth method
DO it=1,itmax
  DO i=dze+1,n*d-die
    IF(check(i)) THEN
      tol=MAX(eps*DCMOD(er(i),ei(i)),eps)
      IF(DABS(ei(i))<tol) THEN
        !update real eigenpair approx
        CALL dgeapprox_eam(p, xr(1,i), xi(1,i), yr(1,i), yi(1,i), er, ei, ncoeff, iseed, berr(i), tol, i, lwork, d, n, check(i))
      ELSE
        !update complex eigenpair approx
        CALL zgeapprox_eam(p, xr(1,i), xi(1,i), yr(1,i), yi(1,i), er, ei, ncoeff, iseed, berr(i), tol, i, lwork, d, n, check(i))
      ENDIF
      IF(.NOT.check(i)) c=c+1
    ENDIF
  ENDDO
  IF(c>=n*d) RETURN
ENDDO
RETURN
END SUBROUTINE dgeeam

!************************************************************************
!			SUBROUTINE DGEAPPROX_EAM			*
!************************************************************************
! Determine whether or not current eigenvalue approximation is close 	*
! enough, if so compute eigenvector, if not update eigenvalue		*
! approximation using Ehrlich-Aberth iterate.	 			*
!************************************************************************
SUBROUTINE dgeapprox_eam(p, xr, xi, yr, yi, er, ei, ncoeff, iseed, berr, tol, i, lwork, d, n, check)
IMPLICIT NONE
!scalar arguments
LOGICAL, INTENT(INOUT) :: check
INTEGER, INTENT(IN) :: d, i, lwork, n
REAL(dp), INTENT(IN) :: tol
REAL(dp), INTENT(INOUT) :: berr
!array arguments
INTEGER, INTENT(INOUT) :: iseed(*)
REAL(dp), INTENT(IN) :: ncoeff(*), p(n,*)
REAL(dp), INTENT(INOUT) :: er(*), ei(*)
REAL(dp), INTENT(INOUT) :: xr(*), xi(*), yr(*), yi(*)
!local scalars
INTEGER :: info
REAL(dp) :: alpha, t
!local arrays
INTEGER, DIMENSION(n) :: jpvt
REAL(dp), DIMENSION(n) :: tau
REAL(dp), DIMENSION(lwork) :: work
REAL(dp), DIMENSION(n,n) :: a
!intrinsic procedures
INTRINSIC :: DABS

!initialize parameters
t=er(i)
jpvt=0
!split into 2 cases
IF(DABS(t)>1) THEN
  t=1/t
  CALL drevgeeval(p, t, a, d, n, 0)
  CALL drevseval(ncoeff, DABS(t), alpha, d, 0)
  a=a/alpha
  CALL dgeqp3(n, n, a, n, jpvt, tau, work, lwork, info)
  IF(DABS(a(n,n))<eps) THEN
    !1st stopping criteria met
    CALL dker1(a, xr, yr, jpvt, tau, work, lwork, n, n, n)
    xi(1:n)=zero; yi(1:n)=zero; ei(i)=zero 
    berr=DABS(a(n,n))
    check=.FALSE.
    RETURN
  ENDIF
  CALL dberrapprox(a, tau, work, iseed, berr, lwork, n)
  IF(berr<eps) THEN
    !2nd stopping criteria met
    CALL dker2(a, xr, yr, jpvt, tau, work, lwork, n)
    xi(1:n)=zero; yi(1:n)=zero
    check=.FALSE.
    RETURN
  ENDIF
  !update eigenvalue approximation
  CALL dgeeamcorr2(p, a, er, ei, tau, work, jpvt, alpha, t, tol, i, lwork, d, n, check)
  IF(.NOT.check) THEN
    !3rd stopping criteria met
    CALL dker2(a, xr, yr, jpvt, tau, work, lwork, n)
    xi(1:n)=zero; yi(1:n)=zero; ei(i)=zero
    RETURN
  ENDIf
ELSE
  CALL dgeeval(p, t, a, d, n, 0)
  CALL dseval(ncoeff, DABS(t), alpha, d, 0)
  a=a/alpha
  CALL dgeqp3(n, n, a, n, jpvt, tau, work, lwork, info)
  IF(DABS(a(n,n))<eps) THEN
    !1st stopping criteria met
    CALL dker1(a, xr, yr, jpvt, tau, work, lwork, n, n, n)
    xi(1:n)=zero; yi(1:n)=zero; ei(i)=zero 
    berr=DABS(a(n,n))
    check=.FALSE.
    RETURN
  ENDIF
  CALL dberrapprox(a, tau, work, iseed, berr, lwork, n)
  IF(berr<eps) THEN
    !2nd stopping criteria met
    CALL dker2(a, xr, yr, jpvt, tau, work, lwork, n)
    xi(1:n)=zero; yi(1:n)=zero; ei(i)=zero
    check=.FALSE.
    RETURN
  ENDIF
  !update eigenvalue approximation
  CALL dgeeamcorr1(p, a, er, ei, tau, work, jpvt, alpha, t, tol, i, lwork, d, n, check)
  IF(.NOT.check) THEN
    !3rd stopping criteria met
    CALL dker2(a, xr, yr, jpvt, tau, work, lwork, n)
    xi(1:n)=zero; yi(1:n)=zero; ei(i)=zero
    RETURN
  ENDIf
ENDIF
RETURN
END SUBROUTINE dgeapprox_eam

!************************************************************************
!			SUBROUTINE DGEEAMCORR1				*
!************************************************************************
! Compute the Laguerre correction term for matrix polynomial p with real*
! coeffs of size n and degree d at real number t, where |t|<=1. All i-1 *
! previously found eigenvalues are stored in (er,ei). The ith eigenvalue*
! approximation is updated in er(i), ei(i).                             *
!************************************************************************
SUBROUTINE dgeeamcorr1(p, a, er, ei, tau, work, jpvt, alpha, t, tol, i, lwork, d, n, check)
IMPLICIT NONE
!scalar arguments
LOGICAL, INTENT(INOUT) :: check
INTEGER, INTENT(IN) :: d, i, lwork, n
REAL(dp), INTENT(IN) :: alpha, tol
REAL(dp), INTENT(INOUT) :: t
!array arguments
INTEGER, INTENT(IN) :: jpvt(*)
REAL(dp), INTENT(IN) :: a(n,n), p(n,*), tau(*)
REAL(dp), INTENT(INOUT) :: er(*), ei(*), work(*)
!local scalars
INTEGER :: info, j, k
REAL(dp) :: temp
COMPLEX(dp) :: x1, y1
!local arrays
REAL(dp), DIMENSION(n) :: v1
REAL(dp), DIMENSION(n,n) :: b
!intrinsic procedures
INTRINSIC :: DBLE, DCMPLX, DIMAG, ZABS

!initiate variables
x1=czero
DO k=1,i-1
  y1=1/DCMPLX(er(i)-er(k),-ei(k))
  x1=x1+y1
ENDDO
DO k=i+1,n*d
  y1=1/DCMPLX(er(i)-er(k),-ei(k))
  x1=x1+y1
ENDDO
!compute b=revp'
CALL dgeeval(p, t, b, d, n, 1)
b=b/alpha
!compute x1=a^{-1}b
CALL dormqr('L','T',n,n,n,a,n,tau,b,n,work,lwork,info)
CALL dtrtrs('U','N','N',n,n,a,n,b,n,info)
DO k=1,n
  v1=b(:,k)
  DO j=1,n
    b(jpvt(j),k)=v1(j)
  ENDDO
ENDDO
!compute y1=p'/p
y1=czero
DO k=1,n
  y1=y1+b(k,k)
ENDDO
!compute Laguerre iterate
x1=y1-x1
y1=1/x1
IF(ZABS(y1)<tol) THEN
  check=.FALSE.
ELSE
  er(i)=er(i)-DBLE(y1)
  ei(i)=-DIMAG(y1)
ENDIF
RETURN
END SUBROUTINE dgeeamcorr1

!************************************************************************
!			SUBROUTINE DGEEAMCORR2				*
!************************************************************************
! Compute the Laguerre correction term for matrix polynomial p with real*
! coeffs of size n and degree d at real number t, where |t|<=1. All i-1 *
! previously found eigenvalues are stored in (er,ei). The ith eigenvalue*
! approximation is updated in er(i), ei(i).                             *
!************************************************************************
SUBROUTINE dgeeamcorr2(p, a, er, ei, tau, work, jpvt, alpha, t, tol, i, lwork, d, n, check)
IMPLICIT NONE
!scalar arguments
LOGICAL, INTENT(INOUT) :: check
INTEGER, INTENT(IN) :: d, i, lwork, n
REAL(dp), INTENT(IN) :: alpha, tol
REAL(dp), INTENT(INOUT) :: t
!array arguments
INTEGER, INTENT(IN) :: jpvt(*)
REAL(dp), INTENT(IN) :: a(n,n), p(n,*), tau(*)
REAL(dp), INTENT(INOUT) :: er(*), ei(*), work(*)
!local scalars
INTEGER :: info, j, k
REAL(dp) :: temp
COMPLEX(dp) :: x1, y1
!local arrays
REAL(dp), DIMENSION(n) :: v1
REAL(dp), DIMENSION(n,n) :: b
!intrinsic procedures
INTRINSIC :: DBLE, DCMPLX, DIMAG, ZABS

!initiate variables
x1=czero
DO k=1,i-1
  y1=1/DCMPLX(er(i)-er(k),-ei(k))
  x1=x1+y1
ENDDO
DO k=i+1,n*d
  y1=1/DCMPLX(er(i)-er(k),-ei(k))
  x1=x1+y1
ENDDO
!compute b=revp'
CALL drevgeeval(p, t, b, d, n, 1)
b=b/alpha
!compute x1=a^{-1}b
CALL dormqr('L','T',n,n,n,a,n,tau,b,n,work,lwork,info)
CALL dtrtrs('U','N','N',n,n,a,n,b,n,info)
DO k=1,n
  v1=b(:,k)
  DO j=1,n
    b(jpvt(j),k)=v1(j)
  ENDDO
ENDDO
!compute y1=p'/p
y1=czero
DO k=1,n
  y1=y1+t*(d-t*b(k,k))
ENDDO
!compute Laguerre iterate
x1=y1-x1
y1=1/x1
IF(ZABS(y1)<tol) THEN
  check=.FALSE.
ELSE
  er(i)=er(i)-DBLE(y1)
  ei(i)=-DIMAG(y1)
ENDIF
RETURN
END SUBROUTINE dgeeamcorr2

!************************************************************************
!			SUBROUTINE ZGEAPPROX_EAM			*
!************************************************************************
! Determine whether or not current eigenvalue approximation is close 	*
! enough, if so compute eigenvector, if not update eigenvalue		*
! approximation using Ehrlich-Aberth iterate.	 			*
!************************************************************************
SUBROUTINE zgeapprox_eam(p, xr, xi, yr, yi, er, ei, ncoeff, iseed, berr, tol, i, lwork, d, n, check)
IMPLICIT NONE
!scalar arguments
LOGICAL, INTENT(INOUT) :: check
INTEGER, INTENT(IN) :: d, i, lwork, n
REAL(dp), INTENT(IN) :: tol
REAL(dp), INTENT(INOUT) :: berr
!array arguments
INTEGER, INTENT(INOUT) :: iseed(*)
REAL(dp), INTENT(IN) :: ncoeff(*), p(n,*)
REAL(dp), INTENT(INOUT) :: er(*), ei(*)
REAL(dp), INTENT(INOUT) :: xr(*), xi(*), yr(*), yi(*)
!local scalars
INTEGER :: info
REAL(dp) :: alpha
COMPLEX(dp) :: t
!local arrays
INTEGER, DIMENSION(n) :: jpvt
REAL(dp), DIMENSION(2*n) :: rwork
COMPLEX(dp), DIMENSION(n) :: tau, x, y
COMPLEX(dp), DIMENSION(lwork) :: work
COMPLEX(dp), DIMENSION(n,n) :: a
!intrinisic procedures
INTRINSIC :: DBLE, DCMPLX, DIMAG, ZABS

!initialize parameters
t=DCMPLX(er(i),ei(i))
jpvt=0
!split into 2 cases
IF(ZABS(t)>1) THEN
  t=1/t
  CALL zrevgeeval(p, t, a, d, n, 0)
  CALL drevseval(ncoeff, ZABS(t), alpha, d, 0)
  a=a/alpha
  CALL zgeqp3(n, n, a, n, jpvt, tau, work, lwork, rwork, info)
  IF(ZABS(a(n,n))<eps) THEN
    !1st stopping criteria met
    CALL zker1(a, x, y, jpvt, tau, work, lwork, n, n, n)
    xr(1:n)=DBLE(x); xi(1:n)=DIMAG(x)
    yr(1:n)=DBLE(y); yi(1:n)=DIMAG(y)
    berr=ZABS(a(n,n))
    check=.FALSE.
    RETURN
  ENDIF
  CALL zberrapprox(a, tau, work, iseed, berr, lwork, n)
  IF(berr<eps) THEN
    !2nd stopping criteria met
    CALL zker2(a, x, y, jpvt, tau, work, lwork, n)
    xr(1:n)=DBLE(x); xi(1:n)=DIMAG(x)
    yr(1:n)=DBLE(y); yi(1:n)=DIMAG(y)
    check=.FALSE.
    RETURN
  ENDIF
  !update eigenvalue approximation
  CALL zgeeamcorr2(p, a, er, ei, tau, work, jpvt, alpha, t, tol, i, lwork, d, n, check)
  IF(.NOT.check) THEN
    !3rd stopping criteria met
    CALL zker2(a, x, y, jpvt, tau, work, lwork, n)
    xr(1:n)=DBLE(x); xi(1:n)=DIMAG(x)
    yr(1:n)=DBLE(y); yi(1:n)=DIMAG(y)
    RETURN
  ENDIf
ELSE
  CALL zgeeval(p, t, a, d, n, 0)
  CALL dseval(ncoeff, ZABS(t), alpha, d, 0)
  a=a/alpha
  CALL zgeqp3(n, n, a, n, jpvt, tau, work, lwork, rwork, info)
  IF(ZABS(a(n,n))<eps) THEN
    !1st stopping criteria met
    CALL zker1(a, x, y, jpvt, tau, work, lwork, n, n, n)
    xr(1:n)=DBLE(x); xi(1:n)=DIMAG(x)
    yr(1:n)=DBLE(y); yi(1:n)=DIMAG(y)
    berr=ZABS(a(n,n))
    check=.FALSE.
    RETURN
  ENDIF
  CALL zberrapprox(a, tau, work, iseed, berr, lwork, n)
  IF(berr<eps) THEN
    !2nd stopping criteria met
    CALL zker2(a, x, y, jpvt, tau, work, lwork, n)
    xr(1:n)=DBLE(x); xi(1:n)=DIMAG(x)
    yr(1:n)=DBLE(y); yi(1:n)=DIMAG(y)
    check=.FALSE.
    RETURN
  ENDIF
  !update eigenvalue approximation
  CALL zgeeamcorr1(p, a, er, ei, tau, work, jpvt, alpha, t, tol, i, lwork, d, n, check)
  IF(.NOT.check) THEN
    !3rd stopping criteria met
    CALL zker2(a, x, y, jpvt, tau, work, lwork, n)
    xr(1:n)=DBLE(x); xi(1:n)=DIMAG(x)
    yr(1:n)=DBLE(y); yi(1:n)=DIMAG(y)
    RETURN
  ENDIf
ENDIF
RETURN
END SUBROUTINE zgeapprox_eam

!************************************************************************
!			SUBROUTINE ZGEEAMCORR1				*
!************************************************************************
! Compute the Laguerre correction term for matrix polynomial p with real*
! coeffs of size n and degree d at real number t, where |t|<=1. All i-1 *
! previously found eigenvalues are stored in (er,ei). The ith eigenvalue*
! approximation is updated in er(i), ei(i).                             *
!************************************************************************
SUBROUTINE zgeeamcorr1(p, a, er, ei, tau, work, jpvt, alpha, t, tol, i, lwork, d, n, check)
IMPLICIT NONE
!scalar arguments
LOGICAL, INTENT(INOUT) :: check
INTEGER, INTENT(IN) :: d, i, lwork, n
REAL(dp), INTENT(IN) :: alpha, tol
COMPLEX(dp), INTENT(INOUT) :: t
!array arguments
INTEGER, INTENT(IN) :: jpvt(*)
REAL(dp), INTENT(IN) :: p(n,*)
REAL(dp), INTENT(INOUT) :: er(*), ei(*)
COMPLEX(dp), INTENT(IN) :: a(n,n), tau(*)
COMPLEX(dp), INTENT(INOUT) :: work(*)
!local scalars
INTEGER :: info, j, k
COMPLEX(dp) :: temp, x1, y1
!local arrays
COMPLEX(dp), DIMENSION(n) :: v1
COMPLEX(dp), DIMENSION(n,n) :: b
!intrinsic procedures
INTRINSIC :: DBLE, DCMPLX, DIMAG, ZABS

!initiate variables
x1=czero
DO k=1,i-1
  y1=1/DCMPLX(er(i)-er(k),ei(i)-ei(k))
  x1=x1+y1
ENDDO
DO k=i+1,n*d
  y1=1/DCMPLX(er(i)-er(k),ei(i)-ei(k))
  x1=x1+y1
ENDDO
!compute b=revp'
CALL zgeeval(p, t, b, d, n, 1)
b=b/alpha
!compute x1=a^{-1}b
CALL zunmqr('L','C',n,n,n,a,n,tau,b,n,work,lwork,info)
CALL ztrtrs('U','N','N',n,n,a,n,b,n,info)
DO k=1,n
  v1=b(:,k)
  DO j=1,n
    b(jpvt(j),k)=v1(j)
  ENDDO
ENDDO
!compute y1=p'/p
y1=czero
DO k=1,n
  y1=y1+b(k,k)
ENDDO
!compute Laguerre iterate
x1=y1-x1
y1=1/x1
IF(ZABS(y1)<tol) THEN
  check=.FALSE.
ELSE
  er(i)=er(i)-DBLE(y1)
  ei(i)=ei(i)-DIMAG(y1)
ENDIF
RETURN
END SUBROUTINE zgeeamcorr1

!************************************************************************
!			SUBROUTINE ZGEEAMCORR2				*
!************************************************************************
! Compute the Laguerre correction term for matrix polynomial p with real*
! coeffs of size n and degree d at real number t, where |t|<=1. All i-1 *
! previously found eigenvalues are stored in (er,ei). The ith eigenvalue*
! approximation is updated in er(i), ei(i).                             *
!************************************************************************
SUBROUTINE zgeeamcorr2(p, a, er, ei, tau, work, jpvt, alpha, t, tol, i, lwork, d, n, check)
IMPLICIT NONE
!scalar arguments
LOGICAL, INTENT(INOUT) :: check
INTEGER, INTENT(IN) :: d, i, lwork, n
REAL(dp), INTENT(IN) :: alpha, tol
COMPLEX(dp), INTENT(INOUT) :: t
!array arguments
INTEGER, INTENT(IN) :: jpvt(*)
REAL(dp), INTENT(IN) :: p(n,*)
REAL(dp), INTENT(INOUT) :: er(*), ei(*)
COMPLEX(dp), INTENT(IN) :: a(n,n), tau(*)
COMPLEX(dp), INTENT(INOUT) :: work(*)
!local scalars
INTEGER :: info, j, k
COMPLEX(dp) :: temp, x1, y1
!local arrays
COMPLEX(dp), DIMENSION(n) :: v1
COMPLEX(dp), DIMENSION(n,n) :: b
!intrinsic procedures
INTRINSIC :: DBLE, DCMPLX, DIMAG, ZABS

!initiate variables
x1=czero
DO k=1,i-1
  y1=1/DCMPLX(er(i)-er(k),ei(i)-ei(k))
  x1=x1+y1
ENDDO
DO k=i+1,n*d
  y1=1/DCMPLX(er(i)-er(k),ei(i)-ei(k))
  x1=x1+y1
ENDDO
!compute b=revp'
CALL zrevgeeval(p, t, b, d, n, 1)
b=b/alpha
!compute x1=a^{-1}b
CALL zunmqr('L','C',n,n,n,a,n,tau,b,n,work,lwork,info)
CALL ztrtrs('U','N','N',n,n,a,n,b,n,info)
DO k=1,n
  v1=b(:,k)
  DO j=1,n
    b(jpvt(j),k)=v1(j)
  ENDDO
ENDDO
!compute y1=p'/p
y1=czero
DO k=1,n
  y1=y1+t*(d-t*b(k,k))
ENDDO
!compute Laguerre iterate
x1=y1-x1
y1=1/x1
IF(ZABS(y1)<tol) THEN
  check=.FALSE.
ELSE
  er(i)=er(i)-DBLE(y1)
  ei(i)=ei(i)-DIMAG(y1)
ENDIF
RETURN
END SUBROUTINE zgeeamcorr2


END MODULE dgeeam_subroutines
