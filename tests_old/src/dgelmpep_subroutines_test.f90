MODULE dgelmpep_subroutines
USE dslmpep_subroutines
IMPLICIT NONE

CONTAINS

!************************************************************************
!			SUBROUTINE DGELMT				*
!************************************************************************
SUBROUTINE dgelmt(p, xr, xi, yr, yi, er, ei, berr, ncoeff, iseed, d, n, opt, ier, iei)
IMPLICIT NONE
!scalar arguments
CHARACTER(LEN=2), INTENT(IN) :: opt
INTEGER, INTENT(IN) :: d, n
!array arguments
INTEGER, INTENT(INOUT) :: iseed(*)
REAL(dp), INTENT(IN) :: ncoeff(*), p(n,*)
REAL(dp), INTENT(INOUT) :: berr(*), er(*), ei(*), ier(*), iei(*)
REAL(dp), INTENT(INOUT) :: xr(n,*), xi(n,*), yr(n,*), yi(n,*)
!local scalars
LOGICAL :: check
INTEGER :: die, dze, i, it, lwork, td
REAL(dp) :: tol
!intrinsic procedures
INTRINSIC :: DABS, MAX

!initial estimates
CALL dgestart(p, xr, xi, yr, yi, er, ei, ncoeff, d, die, dze, lwork, n, opt)
td=n*d-die
ier(1:n*d)=er(1:n*d)
iei(1:n*d)=ei(1:n*d)
!Laguerre's method
DO i=dze+1,td
  DO it=1,itmax
    check=(it==itmax)
    tol=MAX(eps*DCMOD(er(i),ei(i)), eps)
    IF(DABS(ei(i))<tol) THEN
      !update real eigenpair approx
      CALL dgeapprox(p, xr(1,i), xi(1,i), yr(1,i), yi(1,i), er, ei, ncoeff,&
                     iseed, berr(i), tol, i, lwork, d, n, td, check)
    ELSE
      !update complex eigenpair approx
      CALL zgeapprox(p, xr(1,i), xi(1,i), yr(1,i), yi(1,i), er, ei, ncoeff,&
                     iseed, berr(i), tol, i, lwork, d, n, td, check)
    ENDIF
    IF(check) THEN
      EXIT
    ENDIF
  ENDDO
ENDDO
RETURN
END SUBROUTINE dgelmt

!************************************************************************
!			SUBROUTINE DPOSTERRCOND				*
!************************************************************************
! Compute a posteriori backward error, condition, and forward error of	*
! each eigenpair associated with the matrix polynomial p.		*
!************************************************************************
SUBROUTINE dposterrcond(p, xr, xi, yr, yi, er, ei, ncoeff, berr, cond, ferr, d, n)
IMPLICIT NONE
!scalar arguments
INTEGER, INTENT(IN) :: d, n
!array arguments
REAL(dp), INTENT(IN) :: er(*), ei(*), ncoeff(*), p(n,*), xr(n,*), xi(n,*), yr(n,*), yi(n,*)
REAL(dp), INTENT(INOUT) :: berr(*), cond(*), ferr(*)
!local scalars
INTEGER :: i
REAL(dp) :: alpha, bex, bey
COMPLEX(dp) :: t
!local arrays
REAL(dp), DIMENSION(n) :: x
COMPLEX(dp), DIMENSION(n) :: u, v
COMPLEX(dp), DIMENSION(n,n) :: a, b
!external procedures
REAL(dp) :: ddot, dnrm2, dznrm2
COMPLEX(dp) :: zdotc
EXTERNAL :: ddot, dnrm2, dznrm2, zdotc

DO i=1,n*d
  IF(DCMOD(er(i),ei(i))<eps) THEN
    !zero eigenvalues
    CALL dgemv('N', n, n, one, p(1,1), n, xr(1,i), 1, zero, x, 1)
    bex=dnrm2(n, x, 1)/ncoeff(1)
    CALL dgemv('T', n, n, one, p(1,1), n, yr(1,i), 1, zero, x, 1)
    bey=dnrm2(n, x, 1)/ncoeff(1)
    berr(i)=MAX(bex,bey)
    cond(i)=1/DABS(ddot(n, xr(1,i), 1, yr(1,i), 1))
    ferr(i)=berr(i)*cond(i) 
  ELSEIF(DCMOD(er(i),ei(i))<=one) THEN
    !nonzero eigenvalues w/ moduli <=1
    t=DCMPLX(er(i), ei(i))
    CALL zgeeval(p, t, a, d, n, 0)
    CALL dseval(ncoeff, ZABS(t), alpha, d, 0)
    u=DCMPLX(xr(:,i),xi(:,i))
    CALL zgemv('N', n, n, cone, a, n, u, 1, czero, v, 1)
    bex=dznrm2(n, v, 1)
    v=DCMPLX(yr(:,i),yi(:,i))
    CALL zgemv('C', n, n, cone, a, n, v, 1, czero, u, 1)
    bey=dznrm2(n, u, 1)
    berr(i)=MAX(bex,bey)
    CALL zgeeval(p, t, b, d, n, 1)
    u=DCMPLX(xr(:,i),xi(:,i))
    CALL zgemv('N', n, n, cone, b, n, u, 1, czero, v, 1)
    u=DCMPLX(yr(:,i),yi(:,i))
    cond(i)=1/ZABS(t*zdotc(n,u,1,v,1))
    ferr(i)=berr(i)*cond(i)
    berr(i)=berr(i)/alpha
    cond(i)=cond(i)*alpha
  ELSEIF(DCMOD(er(i),ei(i))<big) THEN
    !nonzero eigenvalues w/ moduli >1 and <big
    t=1/DCMPLX(er(i), ei(i))
    CALL zrevgeeval(p, t, a, d, n, 0)
    CALL drevseval(ncoeff, ZABS(t), alpha, d, 0)
    u=DCMPLX(xr(:,i),xi(:,i))
    CALL zgemv('N', n, n, cone, a, n, u, 1, czero, v, 1)
    bex=dznrm2(n, v, 1)
    v=DCMPLX(yr(:,i),yi(:,i))
    CALL zgemv('C', n, n, cone, a, n, v, 1, czero, u, 1)
    bey=dznrm2(n, u, 1)
    berr(i)=MAX(bex,bey)
    CALL zrevgeeval(p, t, b, d, n, 1)
    b=d*a-t*b
    u=DCMPLX(xr(:,i),xi(:,i))
    CALL zgemv('N', n, n, cone, b, n, u, 1, czero, v, 1)
    u=DCMPLX(yr(:,i),yi(:,i))
    cond(i)=1/ZABS(zdotc(n,u,1,v,1))
    ferr(i)=berr(i)*cond(i)
    berr(i)=berr(i)/alpha
    cond(i)=cond(i)*alpha
  ELSE
    !infinite eigenvalues
    CALL dgemv('N', n, n, one, p(1,n*d+1), n, xr(1,i), 1, zero, x, 1)
    bex=dnrm2(n, x, 1)/ncoeff(d+1)
    CALL dgemv('T', n, n, one, p(1,n*d+1), n, yr(1,i), 1, zero, x, 1)
    bey=dnrm2(n, x, 1)/ncoeff(d+1)
    berr(i)=MAX(bex,bey)
    cond(i)=1/DABS(ddot(n, xr(1,i), 1, yr(1,i), 1))
    ferr(i)=berr(i)*cond(i)
  ENDIF
ENDDO
RETURN
END SUBROUTINE dposterrcond

!************************************************************************
!			SUBROUTINE DGELM				*
!************************************************************************
! Compute the eigenvalues and eigenvectors of the matrix polynomial p 	*
! with real coefficients of degree d and size n using Laguerre's method.*
! The norm of the coefficients of p are stored in ncoeff, iseed is used *
! when computing an upper bound on the backward error, stored in berr	*
! of each eigenpair. Eigenvalues are stored in (er,ei) and eigenvectors *
! in (xr,xi) and (yr,yi).						*
!************************************************************************
SUBROUTINE dgelm(p, xr, xi, yr, yi, er, ei, berr, ncoeff, iseed, d, n, opt)
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
LOGICAL :: check
INTEGER :: die, dze, i, it, lwork, td
REAL(dp) :: tol
!intrinsic procedures
INTRINSIC :: DABS, MAX

!initial estimates
CALL dgestart(p, xr, xi, yr, yi, er, ei, ncoeff, d, die, dze, lwork, n, opt)
td=n*d-die
!Laguerre's method
DO i=dze+1,td
  DO it=1,itmax
    check=(it==itmax)
    tol=MAX(eps*DCMOD(er(i),ei(i)), eps)
    IF(DABS(ei(i))<tol) THEN
      !update real eigenpair approx
      CALL dgeapprox(p, xr(1,i), xi(1,i), yr(1,i), yi(1,i), er, ei, ncoeff,&
                     iseed, berr(i), tol, i, lwork, d, n, td, check)
    ELSE
      !update complex eigenpair approx
      CALL zgeapprox(p, xr(1,i), xi(1,i), yr(1,i), yi(1,i), er, ei, ncoeff,&
                     iseed, berr(i), tol, i, lwork, d, n, td, check)
    ENDIF
    IF(check) THEN
      EXIT
    ENDIF
  ENDDO
ENDDO
RETURN
END SUBROUTINE dgelm

!************************************************************************
!			SUBROUTINE DGEAPPROX				*
!************************************************************************
! Determine whether or not current eigenvalue approximation is close 	*
! enough, if so compute eigenvector, if not update eigenvalue		*
! approximation using Laguerre iterate. 				*
!************************************************************************
SUBROUTINE dgeapprox(p, xr, xi, yr, yi, er, ei, ncoeff, iseed, berr, tol, i, lwork, d, n, td, check)
IMPLICIT NONE
!scalar arguments
LOGICAL, INTENT(INOUT) :: check
INTEGER, INTENT(IN) :: d, i, lwork, n, td
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
IF(DABS(t)>one) THEN
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
    check=.TRUE.
    RETURN
  ENDIF
  CALL dberrapprox(a, tau, work, iseed, berr, lwork, n)
  IF(berr<eps) THEN
    !2nd stopping criteria met
    check=.TRUE.
  ENDIF
  IF(check) THEN
    !2nd stopping criteria met, or it==itmax
    CALL dker2(a, xr, yr, jpvt, tau, work, lwork, n)
    xi(1:n)=zero; yi(1:n)=zero; ei(i)=zero
    RETURN
  ENDIF
  !update eigenvalue approximation
  CALL dgelcorr2(p, a, er, ei, tau, work, jpvt, alpha, t, tol, i, lwork, d, n, td, check)
  IF(check) THEN
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
    check=.TRUE.
    RETURN
  ENDIF
  CALL dberrapprox(a, tau, work, iseed, berr, lwork, n)
  IF(berr<eps) THEN
    !2nd stopping criteria met
    check=.TRUE.
  ENDIF
  IF(check) THEN
    !2nd stopping criteria met, or it==itmax
    CALL dker2(a, xr, yr, jpvt, tau, work, lwork, n)
    xi(1:n)=zero; yi(1:n)=zero; ei(i)=zero
    RETURN
  ENDIF
  !update eigenvalue approximation
  CALL dgelcorr1(p, a, er, ei, tau, work, jpvt, alpha, t, tol, i, lwork, d, n, td, check)
  IF(check) THEN
    !3rd stopping criteria met
    CALL dker2(a, xr, yr, jpvt, tau, work, lwork, n)
    xi(1:n)=zero; yi(1:n)=zero; ei(i)=zero
    RETURN
  ENDIf
ENDIF
RETURN
END SUBROUTINE dgeapprox

!************************************************************************
!			SUBROUTINE DGELCORR1				*
!************************************************************************
! Compute the Laguerre correction term for matrix polynomial p with real*
! coeffs of size n and degree d at real number t, where |t|<=1. All i-1 *
! previously found eigenvalues are stored in (er,ei). The ith eigenvalue*
! approximation is updated in er(i), ei(i).                             *
!************************************************************************
SUBROUTINE dgelcorr1(p, a, er, ei, tau, work, jpvt, alpha, t, tol, i, lwork, d, n, td, check)
IMPLICIT NONE
!scalar arguments
LOGICAL, INTENT(INOUT) :: check
INTEGER, INTENT(IN) :: d, i, lwork, n, td
REAL(dp), INTENT(IN) :: alpha, t, tol
!array arguments
INTEGER, INTENT(IN) :: jpvt(*)
REAL(dp), INTENT(IN) :: a(n,n), p(n,*), tau(*)
REAL(dp), INTENT(INOUT) :: er(*), ei(*), work(*)
!local scalars
INTEGER :: info, j, k
REAL(dp) :: temp
COMPLEX(dp) :: x1, x2, y1, y2
!local arrays
REAL(dp), DIMENSION(n) :: v1, v2
REAL(dp), DIMENSION(n,n) :: b, c
!intrinsic procedures
INTRINSIC :: DBLE, DCMPLX, DIMAG, ZABS, ZSQRT

!initiate variables
x1=czero; x2=czero
DO k=1,i-1
  y1=1/DCMPLX(er(i)-er(k),-ei(k))
  x1=x1+y1
  x2=x2+y1**2
ENDDO
!compute b=p', c=p''
CALL dgeeval(p, t, b, d, n, 1)
CALL dgeeval(p, t, c, d, n, 2)
b=b/alpha; c=c/alpha
!compute x1=a^{-1}b and x2=a^{-1}c
CALL dormqr('L','T',n,n,n,a,n,tau,b,n,work,lwork,info)
CALL dormqr('L','T',n,n,n,a,n,tau,c,n,work,lwork,info)
CALL dtrtrs('U','N','N',n,n,a,n,b,n,info)
CALL dtrtrs('U','N','N',n,n,a,n,c,n,info)
DO k=1,n
  v1=b(:,k)
  v2=c(:,k)
  DO j=1,n
    b(jpvt(j),k)=v1(j)
    c(jpvt(j),k)=v2(j)
  ENDDO
ENDDO
!compute y1=p'/p and y2=-(p'/p)'
y1=czero; y2=czero
DO k=1,n
  y1=y1+b(k,k)
  CALL dgemv('N',1,n,one,b(k,1),n,b(1,k),1,zero,temp,1)
  y2=y2+temp-c(k,k)
ENDDO
!compute Laguerre iterate
x1=y1-x1
x2=y2-x2
k=td-i+1
y1=ZSQRT((k-1)*(k*x2-x1**2))
y2=x1-y1; y1=x1+y1
IF(ZABS(y1)>=ZABS(y2)) THEN
  y1=k/y1
  IF(ZABS(y1)<tol) THEN
    !3rd stopping criteria
    check=.TRUE.
  ELSE
    er(i)=er(i)-DBLE(y1)
    ei(i)=-DIMAG(y1)
  ENDIF
ELSE
  y2=k/y2
  IF(ZABS(y2)<tol) THEN
    !3rd stopping criteria
    check=.TRUE.
  ELSE
    er(i)=er(i)-DBLE(y2)
    ei(i)=-DIMAG(y2)
  ENDIF
ENDIF
RETURN
END SUBROUTINE dgelcorr1

!************************************************************************
!			SUBROUTINE DGELCORR2				*
!************************************************************************
! Compute the Laguerre correction term for matrix polynomial p with real*
! coeffs of size n and degree d at real number 1/t, where |t|>1. All i-1*
! previously found eigenvalues are stored in er, ei. The ith eigenvalue *
! approximation is updated in er(i), ei(i).                             *
!************************************************************************
SUBROUTINE dgelcorr2(p, a, er, ei, tau, work, jpvt, alpha, t, tol, i, lwork, d, n, td, check)
IMPLICIT NONE
!scalar arguments
LOGICAL, INTENT(INOUT) :: check
INTEGER, INTENT(IN) :: d, i, lwork, n, td
REAL(dp), INTENT(IN) :: alpha, t, tol
!array arguments
INTEGER, INTENT(IN) :: jpvt(*)
REAL(dp), INTENT(IN) :: a(n,n), p(n,*), tau(*)
REAL(dp), INTENT(INOUT) :: er(*), ei(*), work(*)
!local scalars
INTEGER :: info, j, k
REAL(dp) :: temp
COMPLEX(dp) :: x1, x2, y1, y2
!local arrays
REAL(dp), DIMENSION(n) :: v1, v2
REAL(dp), DIMENSION(n,n) :: b, c
!intrinsic procedures
INTRINSIC :: DBLE, DCMPLX, DIMAG, ZABS, ZSQRT

!initiate variables
x1=czero; x2=czero
DO k=1,i-1
  y1=1/DCMPLX(er(i)-er(k),-ei(k))
  x1=x1+y1
  x2=x2+y1**2
ENDDO
!compute b=revp', c=revp''
CALL drevgeeval(p, t, b, d, n, 1)
CALL drevgeeval(p, t, c, d, n, 2)
b=b/alpha; c=c/alpha
!compute x1=a^{-1}b and x2=a^{-1}c
CALL dormqr('L','T',n,n,n,a,n,tau,b,n,work,lwork,info)
CALL dormqr('L','T',n,n,n,a,n,tau,c,n,work,lwork,info)
CALL dtrtrs('U','N','N',n,n,a,n,b,n,info)
CALL dtrtrs('U','N','N',n,n,a,n,c,n,info)
DO k=1,n
  v1=b(:,k)
  v2=c(:,k)
  DO j=1,n
    b(jpvt(j),k)=v1(j)
    c(jpvt(j),k)=v2(j)
  ENDDO
ENDDO
!compute y1=p'/p and y2=-(p'/p)'
y1=czero; y2=czero
DO k=1,n
  y1=y1+(d-t*b(k,k))
  CALL dgemv('N',1,n,one,b(k,1),n,b(1,k),1,zero,temp,1)
  y2=y2+(d-2*t*b(k,k)+t**2*(temp-c(k,k)))
ENDDO
y1=t*y1
y2=t**2*y2
!compute Laguerre iterate
x1=y1-x1
x2=y2-x2
k=td-i+1
y1=ZSQRT((k-1)*(k*x2-x1**2))
y2=x1-y1; y1=x1+y1
IF(ZABS(y1)>=ZABS(y2)) THEN
  y1=k/y1
  IF(ZABS(y1)<tol) THEN
    !3rd stopping criteria
    check=.TRUE.
  ELSE
    er(i)=er(i)-DBLE(y1)
    ei(i)=-DIMAG(y1)
  ENDIF
ELSE
  y2=k/y2
  IF(ZABS(y2)<tol) THEN
    !3rd stopping criteria
    check=.TRUE.
  ELSE
    er(i)=er(i)-DBLE(y2)
    ei(i)=-DIMAG(y2)
  ENDIF
ENDIF
RETURN
END SUBROUTINE dgelcorr2

!************************************************************************
!			SUBROUTINE DBERRAPPROX				*
!************************************************************************
! Compute an approximation to the backward error associated with the 	*
! smallest eigenvalue of the complex matrix a, which is stored in qr	*
! form (dgeqp3, dgeqrf, dhsqr). Result is stored in berr.		*
!************************************************************************
SUBROUTINE dberrapprox(a, tau, work, iseed, berr, lwork, n)
IMPLICIT NONE
!scalar arguments
INTEGER, INTENT(IN) :: lwork, n
REAL(dp), INTENT(INOUT) :: berr
!array arguments
INTEGER, INTENT(INOUT) :: iseed(*)
REAL(dp), INTENT(IN) :: a(n,n), tau(*)
REAL(dp), INTENT(INOUT) :: work(*)
!local scalars
INTEGER :: info, k
REAL(dp) :: bnorm, xnorm
!local arrays
REAL(dp), DIMENSION(3) :: temp
REAL(dp), DIMENSION(n) :: x
!intrinsic procedures
INTRINSIC :: MINVAL
!external procedures
REAL(dp) :: dnrm2
EXTERNAL :: dnrm2

DO k=1,3
  !random real vector
  CALL dlarnv(2, iseed, n, x)
  bnorm=dnrm2(n, x, 1)
  !solve matrix equation
  CALL dormqr('L','T',n,1,n,a,n,tau,x,n,work,lwork,info)
  CALL dtrtrs('U','N','N',n,1,a,n,x,n,info)
  xnorm=dnrm2(n, x, 1)
  temp(k)=bnorm/xnorm
ENDDO
berr=MINVAL(temp)
RETURN
END SUBROUTINE dberrapprox

!************************************************************************
!			SUBROUTINE ZGEAPPROX				*
!************************************************************************
! Determine whether or not current eigenvalue approximation is close 	*
! enough, if so compute eigenvector, if not update eigenvalue		*
! approximation using Laguerre iterate. 				*
!************************************************************************
SUBROUTINE zgeapprox(p, xr, xi, yr, yi, er, ei, ncoeff, iseed, berr, tol, i, lwork, d, n, td, check)
IMPLICIT NONE
!scalar arguments
LOGICAL, INTENT(INOUT) :: check
INTEGER, INTENT(IN) :: d, i, lwork, n, td
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
t=DCMPLX(er(i), ei(i))
jpvt=0
!split into 2 cases
IF(ZABS(t)>one) THEN
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
    check=.TRUE.
    RETURN
  ENDIF
  CALL zberrapprox(a, tau, work, iseed, berr, lwork, n)
  IF(berr<eps) THEN
    !2nd stopping criteria met
    check=.TRUE.
  ENDIF
  IF(check) THEN
    !2nd stopping criteria met, or it==itmax
    CALL zker2(a, x, y, jpvt, tau, work, lwork, n)
    xr(1:n)=DBLE(x); xi(1:n)=DIMAG(x)
    yr(1:n)=DBLE(y); yi(1:n)=DIMAG(y)
    RETURN
  ENDIF
  !update eigenvalue approximation
  CALL zgelcorr2(p, a, er, ei, tau, work, jpvt, alpha, t, tol, i, lwork, d, n, td, check)
  IF(check) THEN
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
    check=.TRUE.
    RETURN
  ENDIF
  CALL zberrapprox(a, tau, work, iseed, berr, lwork, n)
  IF(berr<eps) THEN
    !2nd stopping criteria met
    check=.TRUE.
  ENDIF
  IF(check) THEN
    !2nd stopping criteria met, or it==itmax
    CALL zker2(a, x, y, jpvt, tau, work, lwork, n)
    xr(1:n)=DBLE(x); xi(1:n)=DIMAG(x)
    yr(1:n)=DBLE(y); yi(1:n)=DIMAG(y)
    RETURN
  ENDIF
  !update eigenvalue approximation
  CALL zgelcorr1(p, a, er, ei, tau, work, jpvt, alpha, t, tol, i, lwork, d, n, td, check)
  IF(check) THEN
    !3rd stopping criteria met
    CALL zker2(a, x, y, jpvt, tau, work, lwork, n)
    xr(1:n)=DBLE(x); xi(1:n)=DIMAG(x)
    yr(1:n)=DBLE(y); yi(1:n)=DIMAG(y)
    RETURN
  ENDIf
ENDIF
RETURN
END SUBROUTINE zgeapprox

!************************************************************************
!			SUBROUTINE ZGELCORR1				*
!************************************************************************
! Compute the Laguerre correction term for matrix polynomial p with real*
! coeffs of size n and degree d at complex number t, where |t|<=1. All 	*
! i-1 previously found eigenvalues are stored in (er,ei). The ith 	*
! eigenvalue approximation is updated in er(i), ei(i).                  *
!************************************************************************
SUBROUTINE zgelcorr1(p, a, er, ei, tau, work, jpvt, alpha, t, tol, i, lwork, d, n, td, check)
IMPLICIT NONE
!scalar arguments
LOGICAL, INTENT(INOUT) :: check
INTEGER, INTENT(IN) :: d, i, lwork, n, td
REAL(dp), INTENT(IN) :: alpha, tol
COMPLEX(dp), INTENT(IN) :: t
!array arguments
INTEGER, INTENT(IN) :: jpvt(*)
REAL(dp), INTENT(IN) :: p(n,*)
REAL(dp), INTENT(INOUT) :: er(*), ei(*)
COMPLEX(dp), INTENT(IN) :: a(n,n), tau(*)
COMPLEX(dp), INTENT(INOUT) :: work(*)
!local scalars
INTEGER :: info, j, k
COMPLEX(dp) :: temp, x1, x2, y1, y2
!local arrays
COMPLEX(dp), DIMENSION(n) :: v1, v2
COMPLEX(dp), DIMENSION(n,n) :: b, c
!intrinsic procedures
INTRINSIC :: DBLE, DCMPLX, DIMAG, ZABS, ZSQRT

!initiate variables
x1=czero; x2=czero
DO k=1,i-1
  y1=1/DCMPLX(er(i)-er(k),ei(i)-ei(k))
  x1=x1+y1
  x2=x2+y1**2
ENDDO
!compute b=revp', c=revp''
CALL zgeeval(p, t, b, d, n, 1)
CALL zgeeval(p, t, c, d, n, 2)
b=b/alpha; c=c/alpha
!compute x1=a^{-1}b and x2=a^{-1}c
CALL zunmqr('L','C',n,n,n,a,n,tau,b,n,work,lwork,info)
CALL zunmqr('L','C',n,n,n,a,n,tau,c,n,work,lwork,info)
CALL ztrtrs('U','N','N',n,n,a,n,b,n,info)
CALL ztrtrs('U','N','N',n,n,a,n,c,n,info)
DO k=1,n
  v1=b(:,k)
  v2=c(:,k)
  DO j=1,n
    b(jpvt(j),k)=v1(j)
    c(jpvt(j),k)=v2(j)
  ENDDO
ENDDO
!compute y1=p'/p and y2=-(p'/p)'
y1=czero; y2=czero
DO k=1,n
  y1=y1+b(k,k)
  CALL zgemv('N',1,n,cone,b(k,1),n,b(1,k),1,czero,temp,1)
  y2=y2+temp-c(k,k)
ENDDO
!compute Laguerre iterate
x1=y1-x1
x2=y2-x2
k=td-i+1
y1=ZSQRT((k-1)*(k*x2-x1**2))
y2=x1-y1; y1=x1+y1
IF(ZABS(y1)>=ZABS(y2)) THEN
  y1=k/y1
  IF(ZABS(y1)<tol) THEN
    !3rd stopping criteria
    check=.TRUE.
  ELSE
    er(i)=er(i)-DBLE(y1)
    ei(i)=ei(i)-DIMAG(y1)
  ENDIF
ELSE
  y2=k/y2
  IF(ZABS(y2)<tol) THEN
    !3rd stopping criteria
    check=.TRUE.
  ELSE
    er(i)=er(i)-DBLE(y2)
    ei(i)=ei(i)-DIMAG(y2)
  ENDIF
ENDIF
RETURN
END SUBROUTINE zgelcorr1

!************************************************************************
!			SUBROUTINE ZGELCORR2				*
!************************************************************************
! Compute the Laguerre correction term for matrix polynomial p with real*
! coeffs of size n and degree d at complex number 1/t, where |t|>1. All *
! i-1 previously found eigenvalues are stored in er, ei. The ith 	*
! eigenvalue approximation is updated in er(i), ei(i).              	*
!************************************************************************
SUBROUTINE zgelcorr2(p, a, er, ei, tau, work, jpvt, alpha, t, tol, i, lwork, d, n, td, check)
IMPLICIT NONE
!scalar arguments
LOGICAL, INTENT(INOUT) :: check
INTEGER, INTENT(IN) :: d, i, lwork, n, td
REAL(dp), INTENT(IN) :: alpha, tol
COMPLEX(dp), INTENT(IN) :: t
!array arguments
INTEGER, INTENT(IN) :: jpvt(*)
REAL(dp), INTENT(IN) :: p(n,*)
REAL(dp), INTENT(INOUT) :: er(*), ei(*)
COMPLEX(dp), INTENT(IN) :: a(n,n), tau(*)
COMPLEX(dp), INTENT(INOUT) :: work(*)
!local scalars
INTEGER :: info, j, k
COMPLEX(dp) :: temp, x1, x2, y1, y2
!local arrays
COMPLEX(dp), DIMENSION(n) :: v1, v2
COMPLEX(dp), DIMENSION(n,n) :: b, c
!intrinsic procedures
INTRINSIC :: DBLE, DCMPLX, DIMAG, ZABS, ZSQRT

!initiate variables
x1=czero; x2=czero
DO k=1,i-1
  y1=1/DCMPLX(er(i)-er(k),ei(i)-ei(k))
  x1=x1+y1
  x2=x2+y1**2
ENDDO
!compute b=revp', c=revp''
CALL zrevgeeval(p, t, b, d, n, 1)
CALL zrevgeeval(p, t, c, d, n, 2)
b=b/alpha; c=c/alpha
!compute x1=a^{-1}b and x2=a^{-1}c
CALL zunmqr('L','C',n,n,n,a,n,tau,b,n,work,lwork,info)
CALL zunmqr('L','C',n,n,n,a,n,tau,c,n,work,lwork,info)
CALL ztrtrs('U','N','N',n,n,a,n,b,n,info)
CALL ztrtrs('U','N','N',n,n,a,n,c,n,info)
DO k=1,n
  v1=b(:,k)
  v2=c(:,k)
  DO j=1,n
    b(jpvt(j),k)=v1(j)
    c(jpvt(j),k)=v2(j)
  ENDDO
ENDDO
!compute y1=p'/p and y2=-(p'/p)'
y1=czero; y2=czero
DO k=1,n
  y1=y1+(d-t*b(k,k))
  CALL zgemv('N',1,n,cone,b(k,1),n,b(1,k),1,czero,temp,1)
  y2=y2+(d-2*t*b(k,k)+t**2*(temp-c(k,k)))
ENDDO
y1=t*y1
y2=t**2*y2
!compute Laguerre iterate
x1=y1-x1
x2=y2-x2
k=td-i+1
y1=ZSQRT((k-1)*(k*x2-x1**2))
y2=x1-y1; y1=x1+y1
IF(ZABS(y1)>=ZABS(y2)) THEN
  y1=k/y1
  IF(ZABS(y1)<tol) THEN
    !3rd stopping criteria
    check=.TRUE.
  ELSE
    er(i)=er(i)-DBLE(y1)
    ei(i)=ei(i)-DIMAG(y1)
  ENDIF
ELSE
  y2=k/y2
  IF(ZABS(y2)<tol) THEN
    !3rd stopping criteria
    check=.TRUE.
  ELSE
    er(i)=er(i)-DBLE(y2)
    ei(i)=ei(i)-DIMAG(y2)
  ENDIF
ENDIF
RETURN
END SUBROUTINE zgelcorr2

!************************************************************************
!			SUBROUTINE ZBERRAPPROX				*
!************************************************************************
! Compute an approximation to the backward error associated with the 	*
! smallest eigenvalue of the complex matrix a, which is stored in qr	*
! form (zgeqp3, zgeqrf, zhsqr). Result is stored in berr.		*
!************************************************************************
SUBROUTINE zberrapprox(a, tau, work, iseed, berr, lwork, n)
IMPLICIT NONE
!scalar arguments
INTEGER, INTENT(IN) :: lwork, n
REAL(dp), INTENT(INOUT) :: berr
!array arguments
INTEGER, INTENT(INOUT) :: iseed(*)
COMPLEX(dp), INTENT(IN) :: a(n,n), tau(*)
COMPLEX(dp), INTENT(INOUT) :: work(*)
!local scalars
INTEGER :: info, k
REAL(dp) :: bnorm, xnorm
!local arrays
REAL(dp), DIMENSION(3) :: temp
COMPLEX(dp), DIMENSION(n) :: x
!intrinsic procedures
INTRINSIC :: MINVAL
!external procedures
REAL(dp) :: dznrm2
EXTERNAL :: dznrm2

DO k=1,3
  !random complex vector
  CALL zlarnv(2, iseed, n, x)
  bnorm=dznrm2(n, x, 1)
  !solve matrix equation
  CALL zunmqr('L','C',n,1,n,a,n,tau,x,n,work,lwork,info)
  CALL ztrtrs('U','N','N',n,1,a,n,x,n,info)
  xnorm=dznrm2(n, x, 1)
  temp(k)=bnorm/xnorm
ENDDO
berr=MINVAL(temp)
RETURN
END SUBROUTINE zberrapprox

!************************************************************************
!			SUBROUTINE DREVGEEVAL				*
!************************************************************************
! Evaluate reversal of scalar polynomial p with real coeffs of degree d,*
! and its der=0,1,2 derivatives at real number 1/t, where |t|>1. 	*
! Returns evaluation in a.						*
!************************************************************************
SUBROUTINE drevgeeval(p, t, a, d, n, der)
IMPLICIT NONE
!scalar arguments
INTEGER, INTENT(IN) :: d, der, n
REAL(dp), INTENT(IN) :: t
!array arguments
REAL(dp), INTENT(IN) :: p(n,*)
REAL(dp), INTENT(INOUT) :: a(n,n)
!local scalars
INTEGER :: k

IF(der==0) THEN
  a=p(:,1:n)
  DO k=2,d+1
    a=t*a+p(:,n*(k-1)+1:n*k)
  ENDDO
ELSEIF(der==1) THEN
  a=d*p(:,1:n)
  DO k=2,d
    a=t*a+(d-k+1)*p(:,n*(k-1)+1:n*k)
  ENDDO
ELSE
  a=d*(d-1)*P(:,1:n)
  DO k=2,d-1
    a=t*a+(d-k+1)*(d-k)*p(:,n*(k-1)+1:n*k)
  ENDDO
ENDIF
RETURN
END SUBROUTINE drevgeeval

!************************************************************************
!			SUBROUTINE DGEEVAL				*
!************************************************************************
! Evaluate general matrix polynomial p of degree d with real coeffs 	*
! of size n, and its der=0,1,2 derivatives at real number t, where	*
! |t|<=1. Returns evaluation in a.					*
!************************************************************************
SUBROUTINE dgeeval(p, t, a, d, n, der)
IMPLICIT NONE
!scalar arguments
INTEGER, INTENT(IN) :: d, der, n
REAL(dp), INTENT(IN) :: t
!array arguments
REAL(dp), INTENT(IN) :: p(n,*)
REAL(dp), INTENT(INOUT) :: a(n,n)
!local scalars
INTEGER :: k

IF(der==0) THEN
  a=p(:,n*d+1:n*(d+1))
  DO k=d,1,-1
    a=t*a+p(:,n*(k-1)+1:n*k)
  ENDDO
ELSEIF(der==1) THEN
  a=d*p(:,n*d+1:n*(d+1))
  DO k=d,2,-1
    a=t*a+(k-1)*p(:,n*(k-1)+1:n*k)
  ENDDO
ELSE
  a=d*(d-1)*p(:,n*d+1:n*(d+1))
  DO k=d,3,-1
    a=t*a+(k-1)*(k-2)*p(:,n*(k-1)+1:n*k)
  ENDDO
ENDIF
RETURN
END SUBROUTINE dgeeval

!************************************************************************
!			SUBROUTINE ZREVGEEVAL				*
!************************************************************************
! Evaluate reversal of scalar polynomial p with real coeffs of degree d,*
! and its der=0,1,2 derivatives at complex number 1/t, where |t|>1. 	*
! Returns evaluation in a.						*
!************************************************************************
SUBROUTINE zrevgeeval(p, t, a, d, n, der)
IMPLICIT NONE
!scalar arguments
INTEGER, INTENT(IN) :: d, der, n
COMPLEX(dp), INTENT(IN) :: t
!array arguments
REAL(dp), INTENT(IN) :: p(n,*)
COMPLEX(dp), INTENT(INOUT) :: a(n,n)
!local scalars
INTEGER :: k

IF(der==0) THEN
  a=p(:,1:n)
  DO k=2,d+1
    a=t*a+p(:,n*(k-1)+1:n*k)
  ENDDO
ELSEIF(der==1) THEN
  a=d*p(:,1:n)
  DO k=2,d
    a=t*a+(d-k+1)*p(:,n*(k-1)+1:n*k)
  ENDDO
ELSE
  a=d*(d-1)*P(:,1:n)
  DO k=2,d-1
    a=t*a+(d-k+1)*(d-k)*p(:,n*(k-1)+1:n*k)
  ENDDO
ENDIF
RETURN
END SUBROUTINE zrevgeeval

!************************************************************************
!			SUBROUTINE ZGEEVAL				*
!************************************************************************
! Evaluate general matrix polynomial p of degree d with real coeffs 	*
! of size n, and its der=0,1,2 derivatives at complex number t, where	*
! |t|<=1. Return evaluation in a.					*
!************************************************************************
SUBROUTINE zgeeval(p, t, a, d, n, der)
IMPLICIT NONE
!scalar arguments
INTEGER, INTENT(IN) :: d, der, n
COMPLEX(dp), INTENT(IN) :: t
!array arguments
REAL(dp), INTENT(IN) :: p(n,*)
COMPLEX(dp), INTENT(INOUT) :: a(n,n)
!local scalars
INTEGER :: k

IF(der==0) THEN
  a=p(:,n*d+1:n*(d+1))
  DO k=d,1,-1
    a=t*a+p(:,n*(k-1)+1:n*k)
  ENDDO
ELSEIF(der==1) THEN
  a=d*p(:,n*d+1:n*(d+1))
  DO k=d,2,-1
    a=t*a+(k-1)*p(:,n*(k-1)+1:n*k)
  ENDDO
ELSE
  a=d*(d-1)*p(:,n*d+1:n*(d+1))
  DO k=d,3,-1
    a=t*a+(k-1)*(k-2)*p(:,n*(k-1)+1:n*k)
  ENDDO
ENDIF
RETURN
END SUBROUTINE zgeeval

!************************************************************************
!			SUBROUTINE DGESTART				*
!************************************************************************
! Compute the initial estimates of the eigenvalues of a matrix          *
! polynomial of degree d with real coeffs of size n. First the dimension*
! of the infinite eigenspace (die) and zero eigenspace (dze) are found. *
! If any, corresponding eigenvectors are stored in (xr,xi) and (yr,yi). *
! Then initial estimates of the remaining finite eigenvalues are found  *
! using either the numerical range (dienr) or the Newton polygon (dienp)*
! and are stored in (er,ei). 						*
!************************************************************************
SUBROUTINE dgestart(p, xr, xi, yr, yi, er, ei, ncoeff, d, die, dze, lwork, n, opt)
IMPLICIT NONE
!scalar arguments
CHARACTER(LEN=2), INTENT(IN) :: opt
INTEGER, INTENT(IN) :: d, n
INTEGER, INTENT(INOUT) :: die, dze, lwork
!array arguments
REAL(dp), INTENT(IN) :: ncoeff(*), p(n,*)
REAL(dp), INTENT(INOUT) :: er(*), ei(*), xr(n,*), xi(n,*), yr(n,*), yi(n,*)
!local scalars
INTEGER :: info, jlo, k
!local arrays
INTEGER, DIMENSION(n) :: jpvt
REAL(dp), DIMENSION(n) :: tau
REAL(dp), DIMENSION(:), ALLOCATABLE :: work
REAL(dp), DIMENSION(n,n) :: a
!intrinsic procedures
INTRINSIC :: DABS, INT

!optimal lwork for dgeqp3 and dormqr
ALLOCATE(work(1))
lwork=-1
a=p(:,n*d+1:n*(d+1))/ncoeff(d+1)
jpvt=0
die=0
CALL dgeqp3(n, n, a, n, jpvt, tau, work, lwork, info)
lwork=INT(work(1))
DEALLOCATE(work)
ALLOCATE(work(lwork))
!infinite eigenspace
CALL dgeqp3(n, n, a, n, jpvt, tau, work, lwork, info)
IF(DABS(a(n,n))<eps) THEN
  jlo=djlo(a,n)
  !update die
  die=n-jlo+1
  !compute eigenvectors
  k=n*d-die+1
  CALL dker(a, xr(1,k), xi(1,k), yr(1,k), yi(1,k), jpvt, tau, work, lwork, jlo, n)
ENDIF
!zero eigenspace
a=p(:,1:n)/ncoeff(1)
jpvt=0
dze=0
CALL dgeqp3(n, n, a, n, jpvt, tau, work, lwork, info)
IF(DABS(a(n,n))<eps) THEN
  jlo=djlo(a,n)
  !update dze
  dze=n-jlo+1
  CALL dker(a, xr(1,1), xi(1,1), yr(1,1), yi(1,1), jpvt, tau, work, lwork, jlo, n)
ENDIF
!initial estimates
IF(opt=='NP') THEN
  CALL dienp(ncoeff, er, ei, d, n)
ELSE
  CALL dienr(p, a, er, ei, tau, work, lwork, d, n)
ENDIF
!store zero and inf eigenvalues
er(1:dze)=zero; ei(1:dze)=zero
er(n*d-die+1:n*d)=big; ei(n*d-die+1:n*d)=zero
!deallocate
DEALLOCATE(work)
RETURN
END SUBROUTINE dgestart

!************************************************************************
!				SUBROUTINE DIENR        		*
!************************************************************************
! Compute initial estimates for the finite eigenvalues of p using the   *
! numerical range of the matrix polynomial. The matrix a and vectors tau*
! and work are used in dormqr to provide orthogonal vectors. The roots 	*
! of the scalar polynomial produced are computed using dslm.		*
!************************************************************************
SUBROUTINE dienr(p, a, er, ei, tau, work, lwork, d, n)
IMPLICIT NONE
!scalar arguments
INTEGER, INTENT(IN) :: d, n
INTEGER, INTENT(INOUT) :: lwork
!array arguments
REAL(dp), INTENT(IN) :: a(n,n), p(n,*), tau(*)
REAL(dp), INTENT(INOUT) :: er(*), ei(*), work(*)
!local scalars
INTEGER :: i, info, j
REAL(dp) :: tr, ti
!local arrays
REAL(dp), DIMENSION(d+1) :: berr, coeff
REAL(dp), DIMENSION(n) :: x, y
!external procedures
REAL(dp) :: ddot
EXTERNAL :: ddot

DO i=1,n
  x=zero; x(i)=one
  CALL dormqr('L', 'N', n, 1, n, a, n, tau, x, n, work, lwork, info)
  DO j=1,d+1
    CALL dgemv('N', n, n, one, p(1,(j-1)*n+1), n, x, 1, zero, y, 1)
    coeff(j)=ddot(n,x,1,y,1)
  ENDDO
  CALL dslm(coeff, er((i-1)*d+1), ei((i-1)*d+1), berr, d)
ENDDO
!sort
DO i=2,n*d
  DO j=i,2,-1
    IF(DCMOD(er(j),ei(j))<DCMOD(er(j-1),ei(j-1))) THEN
      tr=er(j); ti=ei(j)
      er(j)=er(j-1); ei(j)=ei(j-1)
      er(j-1)=tr; ei(j-1)=ti
    ELSE
      EXIT
    ENDIF
  ENDDO
ENDDO
RETURN
END SUBROUTINE dienr

!************************************************************************
!			SUBROUTINE DIENP				*
!************************************************************************
! Compute initial estimates for the finite eigenvalues of p using the   *
! Newton polygon associated with the matrix polynomial.			*
!************************************************************************
SUBROUTINE dienp(ncoeff, er, ei, d, n)
IMPLICIT NONE
!scalar arguments
INTEGER :: d, n
!array arguments
REAL(dp), INTENT(IN) :: ncoeff(*)
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
DO i=1,d+1
  IF(ncoeff(i)>=eps) THEN
    a(i)=DLOG(ncoeff(i))
  ELSE
    a(i)=-one
  ENDIF
ENDDO
!compute upper convex hall
CALL cnvex(d+1,a,h)
!compute initial estimates
iold=1; c=0; th=pi2/d
DO i=2,d+1
  IF(h(i)) THEN
    nzeros=n*(i-iold)
    r=DEXP((a(iold)-a(i))/(i-iold))
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
END SUBROUTINE dienp

!************************************************************************
!				FUNCTION DJLO				*
!************************************************************************
! Find jlo, the smallest integer j such that |a(j,j)|<eps. This function*
! is called only in the case where we know at least one j satisfies.	*
!************************************************************************
FUNCTION djlo(a, n) RESULT(j)
IMPLICIT NONE
!scalar arguments
INTEGER, INTENT(IN) :: n
!array arguments
REAL(dp), INTENT(IN) :: a(n,n)
!local scalars
INTEGER :: j
!intrinsic procedures
INTRINSIC :: DABS

DO j=1,n
  IF(DABS(a(j,j))<eps) THEN
    EXIT
  ENDIF
ENDDO
RETURN
END FUNCTION djlo

!************************************************************************
!				SUBROUTINE DKER				*
!************************************************************************
! Compute left and right kernal vectors associated with real matrix a	*
! of size n. The matrix a comes in column shifted qr form (dgeqp3), with* 
! jpvt, tau, work, lwork. jlo denotes lowest index of a(j,j), such that *
! |a(j,j)|<eps. Results are returned in xr, yr.                         *
!************************************************************************
SUBROUTINE dker(a, xr, xi, yr, yi, jpvt, tau, work, lwork, jlo, n)
IMPLICIT NONE
!scalar arguments
INTEGER, INTENT(IN) :: jlo, lwork, n
!array arguments
INTEGER, INTENT(IN) :: jpvt(*)
REAL(dp), INTENT(IN) :: a(n,n), tau(*)
REAL(dp), INTENT(INOUT) :: work(*), xr(n,*), xi(n,*), yr(n,*), yi(n,*)
!local scalars
INTEGER :: j, k

DO j=jlo,n
  k=1+(n-j)
  CALL dker1(a, xr(1,k), yr(1,k), jpvt, tau, work, lwork, j, jlo, n)
  xi(:,k)=zero; yi(:,k)=zero
ENDDO
RETURN
END SUBROUTINE dker

!************************************************************************
!				SUBROUTINE DKER2			*
!************************************************************************
! Compute left and right eigenvectors associated with the smallest	*
! eigenvalue of the real matrix a of size n. The matrix a comes in	*
! column shifted qr form (dgeqp3), with jpvt, tau, work, lwork. Results *
! are returned in x and y.						*
!************************************************************************
SUBROUTINE dker2(a, x, y, jpvt, tau, work, lwork, n)
IMPLICIT NONE
!scalar arguments
INTEGER, INTENT(IN) :: lwork, n
!array arguments
INTEGER, INTENT(IN) :: jpvt(*)
REAL(dp), INTENT(IN) :: a(n,n), tau(*)
REAL(dp), INTENT(INOUT) :: x(*), y(*), work(*)
!local scalars
INTEGER :: i, info, k
!local arrays
INTEGER, DIMENSION(n) :: jpvt2
REAL(dp), DIMENSION(n) :: temp
!external procedures
REAL(dp) :: dnrm2
EXTERNAL :: dnrm2

!initial estimates
CALL dker1(a, x, y, jpvt, tau, work, lwork, n, n, n)
!jpvt transpose
DO k=1,n
  DO i=1,n
    IF(i==jpvt(k)) THEN
      EXIT
    ENDIF
  ENDDO
  jpvt2(k)=i
ENDDO
!inverse iteration
DO k=1,3
  !===right eigenvector===
  !apply E^{T}
  DO i=1,n
    temp(jpvt2(i))=x(i)
  ENDDO
  !solve (R^{T}R)temp=temp
  CALL dtrtrs('U','T','N',n,1,a,n,temp,n,info)
  CALL dtrtrs('U','N','N',n,1,a,n,temp,n,info)
  !apply E
  DO i=1,n
    x(jpvt(i))=temp(i)
  ENDDO
  !normalize
  x(1:n)=x(1:n)/dnrm2(n,x,1)

  !===left eigenvector===
  !apply Q^{T}
  CALL dormqr('L','T',n,1,n,a,n,tau,y,n,work,lwork,info)
  !solve (RR^{T})y=y
  CALL dtrtrs('U','N','N',n,1,a,n,y,n,info)
  CALL dtrtrs('U','T','N',n,1,a,n,y,n,info)
  !apply Q
  CALL dormqr('L','N',n,1,n,a,n,tau,y,n,work,lwork,info)
  !normalize
  y(1:n)=y(1:n)/dnrm2(n,y,1)
ENDDO
RETURN
END SUBROUTINE dker2

!************************************************************************
!				SUBROUTINE DKER1			*
!************************************************************************
! Compute the jth kernal vector associated with real matrix a of size n.*
! The matrix a comes in column shift qr form (dgeqp3), with jpvt, tau	*
! work, lwork. jlo denotes lowest index of a(j,j), such that		*
! |a(j,j)|<eps, and j denotes current index. Results are returned in x  *
! and y.								*
!************************************************************************
SUBROUTINE dker1(a, x, y, jpvt, tau, work, lwork, j, jlo, n)
IMPLICIT NONE
!scalar arguments
INTEGER, INTENT(IN) :: j, jlo, lwork, n
!array arguments
INTEGER, INTENT(IN) :: jpvt(*)
REAL(dp), INTENT(IN) :: a(n,n), tau(*)
REAL(dp), INTENT(INOUT) :: x(*), y(*), work(*)
!local scalars
INTEGER :: i, info
!external procedures
REAL(dp) :: dnrm2
EXTERNAL :: dnrm2

!===right eigenvector===
x(jlo:n)=zero; x(j)=one
x(1:jlo-1)=-a(1:jlo-1,j)
!solve system ax=x
CALL dtrtrs('U','N','N',jlo-1,1,a,n,x,n,info)
!apply E
DO i=1,n
  y(jpvt(i))=x(i)
ENDDO
x(1:n)=y(1:n)/dnrm2(n,y,1)
!===left eigenvector===
y(1:n)=zero; y(j)=one
!jth column of Q
CALL dormqr('L','N',n,1,n,a,n,tau,y,n,work,lwork,info)
RETURN
END SUBROUTINE dker1

!************************************************************************
!				SUBROUTINE ZKER				*
!************************************************************************
! Compute left and right kernal vectors associated with complex matrix a*
! of size n. The matrix a comes in column shifted qr form (zgeqp3), with* 
! jpvt, tau, work, lwork. jlo denotes lowest index of a(j,j), such that *
! |a(j,j)|<eps. Results are returned in (xr,xi) and (yr,yi).  		*
!************************************************************************
SUBROUTINE zker(a, xr, xi, yr, yi, jpvt, tau, work, lwork, jlo, n)
IMPLICIT NONE
!scalar arguments
INTEGER, INTENT(IN) :: jlo, lwork, n
!array arguments
INTEGER, INTENT(IN) :: jpvt(*)
REAL(dp), INTENT(INOUT) :: xr(n,*), xi(n,*), yr(n,*), yi(n,*)
COMPLEX(dp), INTENT(IN) :: a(n,n), tau(*)
COMPLEX(dp), INTENT(INOUT) :: work(*)
!local scalars
INTEGER :: j, k
!local arrays
COMPLEX(dp), DIMENSION(n) :: x, y
!intrinsic procedures
INTRINSIC :: DBLE, DIMAG

DO j=jlo,n
  k=1+(n-j)
  CALL zker1(a, x, y, jpvt, tau, work, lwork, j, jlo, n)
  !store into xr, xi, yr, and yi
  xr(:,k)=DBLE(x); xi(:,k)=DIMAG(x)
  yr(:,k)=DBLE(y); yi(:,k)=DIMAG(y)
ENDDO
RETURN
END SUBROUTINE zker

!************************************************************************
!				SUBROUTINE ZKER2			*
!************************************************************************
! Compute left and right eigenvectors associated with the smallest	*
! eigenvalue of the complex matrix a of size n. The matrix a comes in	*
! column shifted qr form (zgeqp3), with jpvt, tau, work, lwork. Results *
! are returned in x and y.						*
!************************************************************************
SUBROUTINE zker2(a, x, y, jpvt, tau, work, lwork, n)
IMPLICIT NONE
!scalar arguments
INTEGER, INTENT(IN) :: lwork, n
!array arguments
INTEGER, INTENT(IN) :: jpvt(*)
COMPLEX(dp), INTENT(IN) :: a(n,n), tau(*)
COMPLEX(dp), INTENT(INOUT) :: x(*), y(*), work(*)
!local scalars
INTEGER :: i, info, k
!local arrays
INTEGER, DIMENSION(n) :: jpvt2
COMPLEX(dp), DIMENSION(n) :: temp
!external procedures
REAL(dp) :: dznrm2
EXTERNAL :: dznrm2

!initial estimates
CALL zker1(a, x, y, jpvt, tau, work, lwork, n, n, n)
!PRINT*, 'zker2', ZABS(a(n,n))
!jpvt transpose
DO k=1,n
  DO i=1,n
    IF(i==jpvt(k)) THEN
      EXIT
    ENDIF
  ENDDO
  jpvt2(k)=i
ENDDO
!inverse iteration
DO k=1,3
  !===right eigenvector===
  !apply E^{T}
  DO i=1,n
    temp(jpvt2(i))=x(i)
  ENDDO
  !solve (R*R)x=x
  CALL ztrtrs('U','C','N',n,1,a,n,temp,n,info)
  CALL ztrtrs('U','N','N',n,1,a,n,temp,n,info)
  !apply E
  DO i=1,n
    x(jpvt(i))=temp(i)
  ENDDO
  !normalize
  x(1:n)=x(1:n)/dznrm2(n,x,1)

  !===left eigenvector===
  !apply Q*
  CALL zunmqr('L','C',n,1,n,a,n,tau,y,n,work,lwork,info)
  !solve (RR*)y=y
  CALL ztrtrs('U','N','N',n,1,a,n,y,n,info)
  CALL ztrtrs('U','C','N',n,1,a,n,y,n,info)
  !apply Q
  CALL zunmqr('L','N',n,1,n,a,n,tau,y,n,work,lwork,info)
  !normalize
  y(1:n)=y(1:n)/dznrm2(n,y,1)
ENDDO
RETURN
END SUBROUTINE zker2

!************************************************************************
!				SUBROUTINE ZKER1			*
!************************************************************************
! Compute the jth kernal vector associated with complex matrix a of size*
! n. The matrix a comes in column shift qr form (zgeqp3), with jpvt, tau*
! work, lwork. jlo denotes lowest index of a(j,j), such that		*
! |a(j,j)|<eps, and j denotes current index. Results are returned in x  *
! and y.								*
!************************************************************************
SUBROUTINE zker1(a, x, y, jpvt, tau, work, lwork, j, jlo, n)
!scalar arguments
INTEGER, INTENT(IN) :: j, jlo, lwork, n
!array arguments
INTEGER, INTENT(IN) :: jpvt(*)
COMPLEX(dp), INTENT(IN) :: a(n,n), tau(*)
COMPLEX(dp), INTENT(INOUT) :: x(*), y(*), work(*)
!local scalars
INTEGER :: i, info
!external procedures
REAL(dp) :: dznrm2
EXTERNAL :: dznrm2

!===right eigenvector===
x(jlo:n)=czero; x(j)=cone
x(1:jlo-1)=-a(1:jlo-1,j)
!solve system ax=x
CALL ztrtrs('U','N','N',jlo-1,1,a,n,x,n,info)
!apply E
DO i=1,n
  y(jpvt(i))=x(i)
ENDDO
x(1:n)=y(1:n)/dznrm2(n,y,1)
!===left eigenvector===
y(1:n)=czero; y(j)=cone
!jth column of Q
CALL zunmqr('L','N',n,1,n,a,n,tau,y,n,work,lwork,info)
RETURN
END SUBROUTINE zker1

END MODULE dgelmpep_subroutines
