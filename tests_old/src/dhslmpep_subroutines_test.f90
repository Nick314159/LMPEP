MODULE dhslmpep_subroutines
USE dgelmpep_subroutines
IMPLICIT NONE
REAL(dp), PARAMETER :: two=2.0_dp
COMPLEX(dp), PARAMETER :: ctwo=DCMPLX(two)
CONTAINS

!************************************************************************
!			SUBROUTINE DHSLM				*
!************************************************************************
! Compute the eigenvalues and eigenvectors of the matrix polynomial p 	*
! with real upper Hessenberg coefficients of degree d and size n using 	*
! Laguerre's method. The norm of the coefficients of p are stored in 	*
! ncoeff, iseed is used when computing an upper bound on the backward 	*
! error, stored in berr of each eigenpair. Eigenvalues are stored in 	*
! (er,ei) and eigenvectors in (xr,xi) and (yr,yi).			*
!************************************************************************
SUBROUTINE dhslm(p, xr, xi, yr, yi, er, ei, berr, ncoeff, iseed, d, n)
IMPLICIT NONE
!scalar arguments
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
CALL dgestart(p, xr, xi, yr, yi, er, ei, ncoeff, d, die, dze, lwork, n, 'NR')
td=n*d-die
!Laguerre's method
DO i=dze+1,td
  DO it=1,itmax
    check=(it==itmax)
    tol=MAX(eps*DCMOD(er(i),ei(i)), eps)
    IF(DABS(ei(i))<tol) THEN
      !update real eigenpair approx
      CALL dhsapprox(p, xr(1,i), xi(1,i), yr(1,i), yi(1,i), er, ei, ncoeff,&
                     iseed, berr(i), tol, i, lwork, d, n, td, check)
    ELSE
      !update complex eigenpair approx
      CALL zhsapprox(p, xr(1,i), xi(1,i), yr(1,i), yi(1,i), er, ei, ncoeff,&
                     iseed, berr(i), tol, i, lwork, d, n, td, check)
    ENDIF
    IF(check) THEN
      !PRINT*, i, it, berr(i)
      EXIT
    ENDIF
  ENDDO
ENDDO
RETURN
END SUBROUTINE dhslm

!************************************************************************
!			SUBROUTINE DHSAPPROX				*
!************************************************************************
! Determine whether or not current eigenvalue approximation is close 	*
! enough, if so compute eigenvector, if not update eigenvalue		*
! approximation using Laguerre iterate. 				*
!************************************************************************
SUBROUTINE dhsapprox(p, xr, xi, yr, yi, er, ei, ncoeff, iseed, berr, tol, i, lwork, d, n, td, check)
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
INTEGER :: j, jmax, jmin
REAL(dp) :: alpha, t
!local arrays
REAL(dp), DIMENSION(n-1) :: lda
REAL(dp), DIMENSION(n) :: tau
REAL(dp), DIMENSION(lwork) :: work
REAL(dp), DIMENSION(n,n) :: a
!intrinsic procedures
INTRINSIC :: DABS

!initialize parameters
t=er(i)
!split into 2 cases
IF(DABS(t)>one) THEN
  t=1/t
  CALL drevgeeval(p, t, a, d, n, 0)
  CALL drevseval(ncoeff, DABS(t), alpha, d, 0)
  a=a/alpha
  DO j=1,n-1
    lda(j)=a(j+1,j)
    IF(DABS(lda(j))<eps) THEN
      lda(j)=eps; a(j+1,j)=eps
    ENDIF
  ENDDO
  CALL dhsqr(a, tau, work, n)
  jmax=DHSJMAX(a,n)
  jmin=DHSJMIN(a,n)
  IF(DABS(a(jmin,jmin))<eps) THEN
    !1st stopping criteria met
    CALL dhsker1(a, xr, yr, tau, work, lwork, jmin, jmax, n)
    xi(1:n)=zero; yi(1:n)=zero; ei(i)=zero
    berr=DABS(a(jmin,jmin))
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
    CALL dhsker2(a, xr, yr, tau, work, lwork, jmin, jmax, n)
    xi(1:n)=zero; yi(1:n)=zero; ei(i)=zero
    RETURN
  ENDIF
  !update eigenvalue approximation
  CALL dhslcorr2(p, a, er, ei, lda, tau, work, alpha, t, tol, i, lwork, d, n, td, check)
  IF(check) THEN
    !3rd stopping criteria met
    CALL dhsker2(a, xr, yr, tau, work, lwork, jmin, jmax, n)
    xi(1:n)=zero; yi(1:n)=zero; ei(i)=zero
    RETURN
  ENDIF
ELSE
  CALL dgeeval(p, t, a, d, n, 0)
  CALL dseval(ncoeff, DABS(t), alpha, d, 0)
  a=a/alpha
  DO j=1,n-1
    lda(j)=a(j+1,j)
    IF(DABS(lda(j))<eps) THEN
      lda(j)=eps; a(j+1,j)=eps
    ENDIF
  ENDDO
  CALL dhsqr(a, tau, work, n)
  jmax=DHSJMAX(a,n)
  jmin=DHSJMIN(a,n)
  IF(DABS(a(jmin,jmin))<eps) THEN
    !1st stopping criteria met
    CALL dhsker1(a, xr, yr, tau, work, lwork, jmin, jmax, n)
    xi(1:n)=zero; yi(1:n)=zero; ei(i)=zero
    berr=DABS(a(jmin,jmin))
    check=.TRUE.
    RETURN
  ENDIf
  CALL dberrapprox(a, tau, work, iseed, berr, lwork, n)
  IF(berr<eps) THEN
    !2nd stopping criteria met
    check=.TRUE.
  ENDIF
  IF(check) THEN
    !2nd stopping criteria met, or it==itmax
    CALL dhsker2(a, xr, yr, tau, work, lwork, jmin, jmax, n)
    xi(1:n)=zero; yi(1:n)=zero; ei(i)=zero
    RETURN
  ENDIF
  !update eigenvalue approximation
  CALL dhslcorr1(p, a, er, ei, lda, tau, work, alpha, t, tol, i, lwork, d, n, td, check)
  IF(check) THEN
    !3rd stopping criteria met
    CALL dhsker2(a, xr, yr, tau, work, lwork, jmin, jmax, n)
    xi(1:n)=zero; yi(1:n)=zero; ei(i)=zero
    RETURN
  ENDIF
ENDIF
RETURN
END SUBROUTINE dhsapprox

!************************************************************************
!				SUBROUTINE DHSLCORR1			*
!************************************************************************
! Compute the Laguerre correction term for matrix polynomial p with real*
! upper Hessenberg coeffs of size n and degree d at real number t, where*
! |t|<=1. All i-1 previously found eigenvalues are stored in (er,ei). 	*
! The ith eigenvalue approximation is updated in er(i), ei(i).		*
!************************************************************************
SUBROUTINE dhslcorr1(p, a, er, ei, lda, tau, work, alpha, t, tol, i, lwork, d, n, td, check)
IMPLICIT NONE
!scalar arguments
LOGICAL, INTENT(INOUT) :: check
INTEGER, INTENT(IN) :: d, i, lwork, n, td
REAL(dp), INTENT(IN) :: alpha, t, tol
!array arguments
REAL(dp), INTENT(IN) :: a(n,n), p(n,*), lda(*), tau(*)
REAL(dp), INTENT(INOUT) :: er(*), ei(*), work(*)
!local scalars
INTEGER :: info, k
REAL(dp) :: b0, b1, b2
COMPLEX(dp) :: temp1, temp2, x1, x2, y1, y2
!local arrays
REAL(dp), DIMENSION(n) :: q, v0, v1, v2
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
!compute 1st column of Q^{T} for Hyman's method
q(1)=one; q(2:n)=zero
CALL dormqr('L','T',n,1,n,a,n,tau,q,n,work,lwork,info)
!compute v0 and b0 for Hyman's method
b0=a(n,n)/q(n)
v0(1:n-1)=b0*q(1:n-1)-a(1:n-1,n); v0(n)=one
CALL dtrtrs('U','N','N',n-1,1,a,n,v0,n,info)
!compute v1 and b1 for Hyman's method
CALL dgemv('N',n,n,one,b,n,v0,1,zero,v1,1)
CALL dormqr('L','T',n,1,n,a,n,tau,v1,n,work,lwork,info)
b1=v1(n)/q(n)
v1(1:n-1)=b1*q(1:n-1)-v1(1:n-1); v1(n)=zero
CALL dtrtrs('U','N','N',n-1,1,a,n,v1,n,info)
!compute b2 for Hyman's method
CALL dgemv('N',n,n,one,c,n,v0,1,zero,v2,1)
CALL dgemv('N',n,n,two,b,n,v1,1,one,v2,1)
CALL dormqr('L','T',n,1,n,a,n,tau,v2,n,work,lwork,info)
b2=v2(n)/q(n)

!compute temp1 and temp 2
temp1=czero; temp2=czero
DO k=1,n-1
  y1=b(k+1,k)/lda(k)
  temp1=temp1+y1
  temp2=temp2+c(k+1,k)/lda(k)-y1**2
ENDDO
temp2=temp2+temp1**2
!compute y1=p'/p and y2=-(p'/p)'
y1=b1/b0+temp1
y2=y1**2-(b2/b0+2*(b1/b0)*temp1+temp2)

!deflate
x1=y1-x1
x2=y2-x2
k=td-i+1

!compute Laguerre iterate
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
END SUBROUTINE dhslcorr1

!************************************************************************
!				SUBROUTINE DHSLCORR2			*
!************************************************************************
! Compute the Laguerre correction term for matrix polynomial p with real*
! upper Hessenberg coeffs of size n and degree d at real number 1/t, 	*
! where |t|>1. All i-1 previously found eigenvalues are stored in 	*
! (er,ei). The ith eigenvalue approximation is updated in er(i), ei(i).	*
!************************************************************************
SUBROUTINE dhslcorr2(p, a, er, ei, lda, tau, work, alpha, t, tol, i, lwork, d, n, td, check)
IMPLICIT NONE
!scalar arguments
LOGICAL, INTENT(INOUT) :: check
INTEGER, INTENT(IN) :: d, i, lwork, n, td
REAL(dp), INTENT(IN) :: alpha, t, tol
!array arguments
REAL(dp), INTENT(IN) :: a(n,n), p(n,*), lda(*), tau(*)
REAL(dp), INTENT(INOUT) :: er(*), ei(*), work(*)
!local scalars
INTEGER :: info, k
REAL(dp) :: b0, b1, b2
COMPLEX(dp) :: temp1, temp2, x1, x2, y1, y2
!local arrays
REAL(dp), DIMENSION(n) :: q, v0, v1, v2
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
!compute 1st column of Q^{T} for Hyman's method
q(1)=one; q(2:n)=zero
CALL dormqr('L','T',n,1,n,a,n,tau,q,n,work,lwork,info)
!compute v0 and b0 for Hyman's method
b0=a(n,n)/q(n)
v0(1:n-1)=b0*q(1:n-1)-a(1:n-1,n); v0(n)=one
CALL dtrtrs('U','N','N',n-1,1,a,n,v0,n,info)
!compute v1 and b1 for Hyman's method
CALL dgemv('N',n,n,one,b,n,v0,1,zero,v1,1)
CALL dormqr('L','T',n,1,n,a,n,tau,v1,n,work,lwork,info)
b1=v1(n)/q(n)
v1(1:n-1)=b1*q(1:n-1)-v1(1:n-1); v1(n)=zero
CALL dtrtrs('U','N','N',n-1,1,a,n,v1,n,info)
!compute b2 for Hyman's method
CALL dgemv('N',n,n,one,c,n,v0,1,zero,v2,1)
CALL dgemv('N',n,n,two,b,n,v1,1,one,v2,1)
CALL dormqr('L','T',n,1,n,a,n,tau,v2,n,work,lwork,info)
b2=v2(n)/q(n)

!compute temp1 and temp 2
temp1=czero; temp2=czero
DO k=1,n-1
  y1=b(k+1,k)/lda(k)
  temp1=temp1+y1
  temp2=temp2+c(k+1,k)/lda(k)-y1**2
ENDDO
temp2=temp2+temp1**2
!compute y1=revp'/revp and y2=revp''/revp
y1=b1/b0+temp1
y2=(b2/b0+2*(b1/b0)*temp1+temp2)
!compute y1=p'/p and y2=-(p'/p)'
y2=t**2*(n*d-2*t*y1+t**2*(y1**2-y2))
y1=t*(n*d-t*y1)

!deflate
x1=y1-x1
x2=y2-x2
k=td-i+1

!compute Laguerre iterate
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
END SUBROUTINE dhslcorr2

!************************************************************************
!				SUBROUTINE DHSQR			*
!************************************************************************
! Compute QR factorization of real upper Hessenberg matrix a. Tau stores*
! elementary reflectors and upper triangulr result is tored in a.	*
!************************************************************************
SUBROUTINE dhsqr(a, tau, work, n)
IMPLICIT NONE
!scalar arguments
INTEGER, INTENT(IN) :: n
!array arguments
REAL(dp), INTENT(INOUT) :: a(n,n), tau(*), work(*)
!local scalars
INTEGER :: j
REAL(dp) :: ajj

DO j=1,n-1
  !generate elementary reflector
  CALL dlarfg(2, a(j,j), a(j+1,j), 1, tau(j))
  !apply elementary reflector
  ajj=a(j,j)
  a(j,j)=one
  CALL dlarf('L', 2, n-j, a(j,j), 1, tau(j), a(j,j+1), n, work)
  a(j,j)=ajj
ENDDO
!nth reflector
CALL dlarfg(1, a(n,n), a(n,n), 1, tau(n))
RETURN
END SUBROUTINE dhsqr

!************************************************************************
!				FUNCTION DHSJMIN			*
!************************************************************************
! Compute j such that |a(j,j)| is minimized, or first index such that	* 
! |a(j,j)|<eps.				 				*
!************************************************************************
FUNCTION dhsjmin(a, n) RESULT(jmin)
IMPLICIT NONE
!scalar arguments
INTEGER, INTENT(IN) :: n
!array arguments
REAL(dp), INTENT(IN) :: a(n,n)
!local scalars
INTEGER :: j, jmin
!intrinsic procedures
INTRINSIC :: DABS

jmin=1
DO j=2,n
  IF(DABS(a(jmin,jmin))<eps) RETURN

  IF(DABS(a(j,j))<DABS(a(jmin,jmin))) THEN
    jmin=j
  ENDIF
ENDDO
RETURN
END FUNCTION dhsjmin

!************************************************************************
!				FUNCTION DHSJMAX			*
!************************************************************************
! Compute j such that |a(j,j)| is minimized, or last index such that	* 
! |a(j,j)|<eps.				 				*
!************************************************************************
FUNCTION dhsjmax(a, n) RESULT(jmax)
IMPLICIT NONE
!scalar arguments
INTEGER, INTENT(IN) :: n
!array arguments
REAL(dp), INTENT(IN) :: a(n,n)
!local scalars
INTEGER :: j, jmax
!intrinsic procedures
INTRINSIC :: DABS

jmax=n
DO j=n-1,1,-1
  IF(DABS(a(jmax,jmax))<eps) RETURN

  IF(DABS(a(j,j))<DABS(a(jmax,jmax))) THEN
    jmax=j
  ENDIF
ENDDO
RETURN
END FUNCTION dhsjmax

!************************************************************************
!			SUBROUTINE ZHSAPPROX				*
!************************************************************************
! Determine whether or not current eigenvalue approximation is close 	*
! enough, if so compute eigenvector, if not update eigenvalue		*
! approximation using Laguerre iterate. 				*
!************************************************************************
SUBROUTINE zhsapprox(p, xr, xi, yr, yi, er, ei, ncoeff, iseed, berr, tol, i, lwork, d, n, td, check)
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
INTEGER :: j, jmax, jmin
REAL(dp) :: alpha
COMPLEX(dp) :: t
!local arrays
COMPLEX(dp), DIMENSION(n-1) :: lda
COMPLEX(dp), DIMENSION(n) :: tau, x, y
COMPLEX(dp), DIMENSION(lwork) :: work
COMPLEX(dp), DIMENSION(n,n) :: a
!intrinsic procedures
INTRINSIC :: DBLE, DCMPLX, DIMAG, ZABS

!initialize parameters
t=DCMPLX(er(i),ei(i))
!split into 2 cases
IF(ZABS(t)>one) THEN
  t=1/t
  CALL zrevgeeval(p, t, a, d, n, 0)
  CALL drevseval(ncoeff, ZABS(t), alpha, d, 0)
  a=a/alpha
  DO j=1,n-1
    lda(j)=a(j+1,j)
    IF(ZABS(lda(j))<eps) THEN
      lda(j)=eps; a(j+1,j)=eps
    ENDIF
  ENDDO
  CALL zhsqr(a, tau, work, n)
  jmax=ZHSJMAX(a,n)
  jmin=ZHSJMIN(a,n)
  IF(ZABS(a(jmin,jmin))<eps) THEN
    !1st stopping criteria met
    CALL zhsker1(a, x, y, tau, work, lwork, jmin, jmax, n)
    xr(1:n)=DBLE(x); xi(1:n)=DIMAG(x)
    yr(1:n)=DBLE(y); yi(1:n)=DIMAG(y)
    berr=ZABS(a(jmin,jmin))
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
    CALL zhsker2(a, x, y, tau, work, lwork, jmin, jmax, n)
    xr(1:n)=DBLE(x); xi(1:n)=DIMAG(x)
    yr(1:n)=DBLE(y); yi(1:n)=DIMAG(y)
    RETURN
  ENDIF
  !update eigenvalue approximation
  CALL zhslcorr2(p, a, er, ei, lda, tau, work, alpha, t, tol, i, lwork, d, n, td, check)
  IF(check) THEN
    !3rd stopping criteria met
    CALL zhsker2(a, x, y, tau, work, lwork, jmin, jmax, n)
    xr(1:n)=DBLE(x); xi(1:n)=DIMAG(x)
    yr(1:n)=DBLE(y); yi(1:n)=DIMAG(y)
    RETURN
  ENDIF
ELSE
  CALL zgeeval(p, t, a, d, n, 0)
  CALL dseval(ncoeff, ZABS(t), alpha, d, 0)
  a=a/alpha
  DO j=1,n-1
    lda(j)=a(j+1,j)
    IF(ZABS(lda(j))<eps) THEN
      lda(j)=eps; a(j+1,j)=eps
    ENDIF
  ENDDO
  CALL zhsqr(a, tau, work, n)
  jmax=ZHSJMAX(a,n)
  jmin=ZHSJMIN(a,n)
  IF(ZABS(a(jmin,jmin))<eps) THEN
    !1st stopping criteria met
    CALL zhsker1(a, x, y, tau, work, lwork, jmin, jmax, n)
    xr(1:n)=DBLE(x); xi(1:n)=DIMAG(x)
    yr(1:n)=DBLE(y); yi(1:n)=DIMAG(y)
    berr=ZABS(a(jmin,jmin))
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
    CALL zhsker2(a, x, y, tau, work, lwork, jmin, jmax, n)
    xr(1:n)=DBLE(x); xi(1:n)=DIMAG(x)
    yr(1:n)=DBLE(y); yi(1:n)=DIMAG(y)
    RETURN
  ENDIF
  !update eigenvalue approximation
  CALL zhslcorr1(p, a, er, ei, lda, tau, work, alpha, t, tol, i, lwork, d, n, td, check)
  IF(check) THEN
    !3rd stopping criteria met
    CALL zhsker2(a, x, y, tau, work, lwork, jmin, jmax, n)
    xr(1:n)=DBLE(x); xi(1:n)=DIMAG(x)
    yr(1:n)=DBLE(y); yi(1:n)=DIMAG(y)
    RETURN
  ENDIF
ENDIF
RETURN
END SUBROUTINE zhsapprox

!************************************************************************
!				SUBROUTINE ZHSLCORR1			*
!************************************************************************
! Compute the Laguerre correction term for matrix polynomial p with	*
! complex upper Hessenberg coeffs of size n and degree d at complex 	*
! number t, where |t|<=1. All i-1 previously found eigenvalues are 	*
! stored in (er,ei). The ith eigenvalue approximation is updated in 	*
! er(i), ei(i).								*
!************************************************************************
SUBROUTINE zhslcorr1(p, a, er, ei, lda, tau, work, alpha, t, tol, i, lwork, d, n, td, check)
IMPLICIT NONE
!scalar arguments
LOGICAL, INTENT(INOUT) :: check
INTEGER, INTENT(IN) :: d, i, lwork, n, td
REAL(dp), INTENT(IN) :: alpha, tol
COMPLEX(dp), INTENT(IN) :: t
!array arguments
REAL(dp), INTENT(IN) :: p(n,*)
REAL(dp), INTENT(INOUT) :: er(*), ei(*)
COMPLEX(dp), INTENT(IN) :: a(n,n), lda(*), tau(*)
COMPLEX(dp), INTENT(INOUT) :: work(*)
!local scalars
INTEGER :: info, k
COMPLEX(dp) :: b0, b1, b2, temp1, temp2, x1, x2, y1, y2
!local arrays
COMPLEX(dp), DIMENSION(n) :: q, v0, v1, v2
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
!compute b=p', c=p''
CALL zgeeval(p, t, b, d, n, 1)
CALL zgeeval(p, t, c, d, n, 2)
b=b/alpha; c=c/alpha
!compute 1st column of Q^{*} for Hyman's method
q(1)=cone; q(2:n)=czero
CALL zunmqr('L','C',n,1,n,a,n,tau,q,n,work,lwork,info)
!compute v0 and b0 for Hyman's method
b0=a(n,n)/q(n)
v0(1:n-1)=b0*q(1:n-1)-a(1:n-1,n); v0(n)=cone
CALL ztrtrs('U','N','N',n-1,1,a,n,v0,n,info)
!compute v1 and b1 for Hyman's method
CALL zgemv('N',n,n,cone,b,n,v0,1,czero,v1,1)
CALL zunmqr('L','C',n,1,n,a,n,tau,v1,n,work,lwork,info)
b1=v1(n)/q(n)
v1(1:n-1)=b1*q(1:n-1)-v1(1:n-1); v1(n)=czero
CALL ztrtrs('U','N','N',n-1,1,a,n,v1,n,info)
!compute b2 for Hyman's method
CALL zgemv('N',n,n,cone,c,n,v0,1,czero,v2,1)
CALL zgemv('N',n,n,ctwo,b,n,v1,1,cone,v2,1)
CALL zunmqr('L','C',n,1,n,a,n,tau,v2,n,work,lwork,info)
b2=v2(n)/q(n)

!compute temp1 and temp 2
temp1=czero; temp2=czero
DO k=1,n-1
  y1=b(k+1,k)/lda(k)
  temp1=temp1+y1
  temp2=temp2+c(k+1,k)/lda(k)-y1**2
ENDDO
temp2=temp2+temp1**2
!compute y1=p'/p and y2=-(p'/p)'
y1=b1/b0+temp1
y2=y1**2-(b2/b0+2*(b1/b0)*temp1+temp2)

!deflate
x1=y1-x1
x2=y2-x2
k=td-i+1

!compute Laguerre iterate
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
END SUBROUTINE zhslcorr1

!************************************************************************
!				SUBROUTINE ZHSLCORR2			*
!************************************************************************
! Compute the Laguerre correction term for matrix polynomial p with	*
! complex upper Hessenberg coeffs of size n and degree d at complex 	*
! number 1/t, where |t|>1. All i-1 previously found eigenvalues are 	*
! stored in (er,ei). The ith eigenvalue approximation is updated in 	*
! er(i), ei(i).								*
!************************************************************************
SUBROUTINE zhslcorr2(p, a, er, ei, lda, tau, work, alpha, t, tol, i, lwork, d, n, td, check)
IMPLICIT NONE
!scalar arguments
LOGICAL, INTENT(INOUT) :: check
INTEGER, INTENT(IN) :: d, i, lwork, n, td
REAL(dp), INTENT(IN) :: alpha, tol
COMPLEX(dp), INTENT(IN) :: t
!array arguments
REAL(dp), INTENT(IN) :: p(n,*)
REAL(dp), INTENT(INOUT) :: er(*), ei(*)
COMPLEX(dp), INTENT(IN) :: a(n,n), lda(*), tau(*)
COMPLEX(dp), INTENT(INOUT) :: work(*)
!local scalars
INTEGER :: info, k
COMPLEX(dp) :: b0, b1, b2, temp1, temp2, x1, x2, y1, y2
!local arrays
COMPLEX(dp), DIMENSION(n) :: q, v0, v1, v2
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
!compute 1st column of Q^{*} for Hyman's method
q(1)=cone; q(2:n)=czero
CALL zunmqr('L','C',n,1,n,a,n,tau,q,n,work,lwork,info)
!compute v0 and b0 for Hyman's method
b0=a(n,n)/q(n)
v0(1:n-1)=b0*q(1:n-1)-a(1:n-1,n); v0(n)=cone
CALL ztrtrs('U','N','N',n-1,1,a,n,v0,n,info)
!compute v1 and b1 for Hyman's method
CALL zgemv('N',n,n,cone,b,n,v0,1,czero,v1,1)
CALL zunmqr('L','C',n,1,n,a,n,tau,v1,n,work,lwork,info)
b1=v1(n)/q(n)
v1(1:n-1)=b1*q(1:n-1)-v1(1:n-1); v1(n)=czero
CALL ztrtrs('U','N','N',n-1,1,a,n,v1,n,info)
!compute b2 for Hyman's method
CALL zgemv('N',n,n,cone,c,n,v0,1,czero,v2,1)
CALL zgemv('N',n,n,ctwo,b,n,v1,1,cone,v2,1)
CALL zunmqr('L','C',n,1,n,a,n,tau,v2,n,work,lwork,info)
b2=v2(n)/q(n)

!compute temp1 and temp 2
temp1=czero; temp2=czero
DO k=1,n-1
  y1=b(k+1,k)/lda(k)
  temp1=temp1+y1
  temp2=temp2+c(k+1,k)/lda(k)-y1**2
ENDDO
temp2=temp2+temp1**2
!compute y1=revp'/revp and y2=revp''/revp
y1=b1/b0+temp1
y2=(b2/b0+2*(b1/b0)*temp1+temp2)
!compute y1=p'/p and y2=-(p'/p)'
y2=t**2*(n*d-2*t*y1+t**2*(y1**2-y2))
y1=t*(n*d-t*y1)

!deflate
x1=y1-x1
x2=y2-x2
k=td-i+1

!compute Laguerre iterate
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
END SUBROUTINE zhslcorr2

!************************************************************************
!				SUBROUTINE ZHSQR			*
!************************************************************************
! Compute QR factorization of complex upper Hessenberg matrix a. Tau 	*
! stores elementary reflectors and upper triangulr result is tored in a.*
!************************************************************************
SUBROUTINE zhsqr(a, tau, work, n)
IMPLICIT NONE
!scalar arguments
INTEGER, INTENT(IN) :: n
!array arguments
COMPLEX(dp), INTENT(INOUT) :: a(n,n), tau(*), work(*)
!local scalars
INTEGER :: j
COMPLEX(dp) :: ajj

DO j=1,n-1
  !generate elementary reflector
  CALL zlarfg(2, a(j,j), a(j+1,j), 1, tau(j))
  !apply elementary reflector
  ajj=a(j,j)
  a(j,j)=cone
  CALL zlarf('L', 2, n-j, a(j,j), 1, DCONJG(tau(j)), a(j,j+1), n, work)
  a(j,j)=ajj
ENDDO
!nth reflector
CALL zlarfg(1, a(n,n), a(n,n), 1, tau(n))
RETURN
END SUBROUTINE zhsqr

!************************************************************************
!				FUNCTION ZHSJMIN			*
!************************************************************************
! Compute j such that |a(j,j)| is minimized, or first index such that	* 
! |a(j,j)|<eps.				 				*
!************************************************************************
FUNCTION zhsjmin(a, n) RESULT(jmin)
IMPLICIT NONE
!scalar arguments
INTEGER, INTENT(IN) :: n
!array arguments
COMPLEX(dp), INTENT(IN) :: a(n,n)
!local scalars
INTEGER :: j, jmin
!intrinsic procedures
INTRINSIC :: ZABS

jmin=1
DO j=2,n
  IF(ZABS(a(jmin,jmin))<eps) RETURN

  IF(ZABS(a(j,j))<ZABS(a(jmin,jmin))) THEN
    jmin=j
  ENDIF
ENDDO
RETURN
END FUNCTION zhsjmin

!************************************************************************
!				FUNCTION ZHSJMAX			*
!************************************************************************
! Compute j such that |a(j,j)| is minimized, or last index such that	* 
! |a(j,j)|<eps.				 				*
!************************************************************************
FUNCTION zhsjmax(a, n) RESULT(jmax)
IMPLICIT NONE
!scalar arguments
INTEGER, INTENT(IN) :: n
!array arguments
COMPLEX(dp), INTENT(IN) :: a(n,n)
!local scalars
INTEGER :: j, jmax
!intrinsic procedures
INTRINSIC :: ZABS

jmax=n
DO j=n-1,1,-1
  IF(ZABS(a(jmax,jmax))<eps) RETURN

  IF(ZABS(a(j,j))<ZABS(a(jmax,jmax))) THEN
    jmax=j
  ENDIF
ENDDO
RETURN
END FUNCTION zhsjmax

!************************************************************************
!				SUBROUTINE DHSKER2			*
!************************************************************************
! Compute left and right singular vectors associated with the smallest	*
! eigenvalue of the real upper Hessenberg matrix a of size n. The matrix*
! a comes in qr form (dhsqr) with tau, work, and lwork. Results are	*
! returned in x and y.							*
!************************************************************************
SUBROUTINE dhsker2(a, x, y, tau, work, lwork, jmin, jmax, n)
IMPLICIT NONE
!scalar arguments
INTEGER, INTENT(IN) :: jmax, jmin, lwork, n
!array arguments
REAL(dp), INTENT(IN) :: a(n,n), tau(*)
REAL(dp), INTENT(INOUT) :: x(*), y(*), work(*)
!local scalars
INTEGER :: info, k
!external procedures
REAL(dp) :: dnrm2
EXTERNAL :: dnrm2

!initial estimates
CALL dhsker1(a, x, y, tau, work, lwork, jmin, jmax, n)
!inverse iteration
DO k=1,3
  !===right eigenvector===
  !solve (R^{T}R)x=x
  CALL dtrtrs('U','T','N',n,1,a,n,x,n,info)
  CALL dtrtrs('U','N','N',n,1,a,n,x,n,info)
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
END SUBROUTINE dhsker2

!************************************************************************
!				SUBROUTINE DHSKER1			*
!************************************************************************
! Computer the kernal vector associated with real upper Hessenberg	*
! matrix of size n. The matrix a comes in qr form (dhsqr), with tau and	*
! work. jmin denotes the index returned by DHSJMIN.			*
!************************************************************************
SUBROUTINE dhsker1(a, x, y, tau, work, lwork, jmin, jmax, n)
IMPLICIT NONE
!scalar arguments
INTEGER, INTENT(IN) :: jmax, jmin, lwork, n
!array arguments
REAL(dp), INTENT(IN) :: a(n,n), tau(*)
REAL(dp), INTENT(INOUT) :: x(*), y(*), work(*)
!local scalars
INTEGER :: info
!external procedures
REAL(dp) :: dnrm2
EXTERNAL :: dnrm2

!===right eigenvector===
IF(jmin==1) THEN
  x(1)=one; x(2:n)=zero
ELSE
  x(jmin)=one; x(jmin+1:n)=zero
  x(1:jmin-1)=-a(1:jmin-1,jmin)
  !solve system ax=x
  CALL dtrtrs('U','N','N',jmin-1,1,a,n,x,n,info)
  x(1:n)=x(1:n)/dnrm2(n,x,1)
ENDIF
!===left eigenvector===
IF(jmax==n) THEN
  y(1:n-1)=zero; y(n)=one
  CALL dormqr('L','N',n,1,n,a,n,tau,y,n,work,lwork,info)
ELSE
  y(1:jmax-1)=zero; y(jmax)=one
  y(jmax+1:n)=-a(jmax,jmax+1:n)
  !solve system R^{T}y=y
  CALL dtrtrs('U','T','N',n-jmax,1,a(jmax+1,jmax+1),n,y(jmax+1),n,info)
  !apply Q
  CALL dormqr('L','N',n,1,n,a,n,tau,y,n,work,lwork,info)
  y(1:n)=y(1:n)/dnrm2(n,y,1)
ENDIF
RETURN
END SUBROUTINE dhsker1

!************************************************************************
!				SUBROUTINE ZHSKER2			*
!************************************************************************
! Compute left and right singular vectors associated with the smallest	*
! eigenvalue of the real upper Hessenberg matrix a of size n. The matrix*
! a comes in qr form (dhsqr) with tau, work, and lwork. Results are	*
! returned in x and y.							*
!************************************************************************
SUBROUTINE zhsker2(a, x, y, tau, work, lwork, jmin, jmax, n)
IMPLICIT NONE
!scalar arguments
INTEGER, INTENT(IN) :: jmax, jmin, lwork, n
!array arguments
COMPLEX(dp), INTENT(IN) :: a(n,n), tau(*)
COMPLEX(dp), INTENT(INOUT) :: x(*), y(*), work(*)
!local scalars
INTEGER :: info, k
!external procedures
REAL(dp) :: dznrm2
EXTERNAL :: dznrm2

!initial estimates
CALL zhsker1(a, x, y, tau, work, lwork, jmin, jmax, n)
!inverse iteration
DO k=1,3
  !===right eigenvector===
  !solve (R^{*}R)x=x
  CALL ztrtrs('U','C','N',n,1,a,n,x,n,info)
  CALL ztrtrs('U','N','N',n,1,a,n,x,n,info)
  !normalize
  x(1:n)=x(1:n)/dznrm2(n,x,1)
  
  !===left eigenvector===
  !apply Q^{*}
  CALL zunmqr('L','C',n,1,n,a,n,tau,y,n,work,lwork,info)
  !solve (RR^{*})y=y
  CALL ztrtrs('U','N','N',n,1,a,n,y,n,info)
  CALL ztrtrs('U','C','N',n,1,a,n,y,n,info)
  !apply Q
  CALL zunmqr('L','N',n,1,n,a,n,tau,y,n,work,lwork,info)
  !normalize
  y(1:n)=y(1:n)/dznrm2(n,y,1)
ENDDO
RETURN
END SUBROUTINE zhsker2

!************************************************************************
!				SUBROUTINE ZHSKER1			*
!************************************************************************
! Computer the kernal vector associated with complex upper Hessenberg	*
! matrix of size n. The matrix a comes in qr form (zhsqr), with tau and	*
! work. jmin denotes the index returned by DHSJMIN.			*
!************************************************************************
SUBROUTINE zhsker1(a, x, y, tau, work, lwork, jmin, jmax, n)
IMPLICIT NONE
!scalar arguments
INTEGER, INTENT(IN) :: jmax, jmin, lwork, n
!array arguments
COMPLEX(dp), INTENT(IN) :: a(n,n), tau(*)
COMPLEX(dp), INTENT(INOUT) :: x(*), y(*), work(*)
!local scalars
INTEGER :: info
!external procedures
REAL(dp) :: dznrm2
EXTERNAL :: dznrm2

!===right eigenvector===
IF(jmin==1) THEN
  x(1)=cone; x(2:n)=czero
ELSE
  x(jmin)=cone; x(jmin+1:n)=czero
  x(1:jmin-1)=-a(1:jmin-1,jmin)
  !solve system ax=x
  CALL ztrtrs('U','N','N',jmin-1,1,a,n,x,n,info)
  x(1:n)=x(1:n)/dznrm2(n,x,1)
ENDIF
!===left eigenvector===
IF(jmax==n) THEN
  y(1:n-1)=czero; y(n)=cone
  !apply Q
  CALL zunmqr('L','N',n,1,n,a,n,tau,y,n,work,lwork,info)
ELSE
  y(1:jmax-1)=czero; y(jmax)=cone
  y(jmax+1:n)=-DCONJG(a(jmax,jmax+1:n))
  !solve system R*y=y
  CALL ztrtrs('U','C','N',n-jmax,1,a(jmax+1,jmax+1),n,y(jmax+1),n,info)
  !apply Q
  CALL zunmqr('L','N',n,1,n,a,n,tau,y,n,work,lwork,info)
  y(1:n)=y(1:n)/dznrm2(n,y,1)
ENDIF
RETURN
END SUBROUTINE zhsker1

END MODULE dhslmpep_subroutines
