PROGRAM dgelmpep_driver
IMPLICIT NONE
INTEGER, PARAMETER :: dp=KIND(0.0D0)
REAL(dp), PARAMETER :: eps=EPSILON(0.0_dp), big=HUGE(0.0_dp)
REAL(dp), PARAMETER :: zero=0.0_dp, one=1.0_dp
COMPLEX(dp), PARAMETER :: czero=DCMPLX(zero), cone=DCMPLX(one)

!solve problem using QUADEIG
CHARACTER(LEN=1) :: jobVL, jobVR
INTEGER :: d, n, opt, info, iwarn
REAL(dp) :: tol
REAL(dp), DIMENSION(:), ALLOCATABLE :: x, work
REAL(dp), DIMENSION(:), ALLOCATABLE :: alphar, alphai, beta, s, beVl, beVR, er, ei, berr, cond, ferr, ncoeff
REAL(dp), DIMENSION(:,:), ALLOCATABLE :: A, B, C, VL, VR, p, xr, xi, yr, yi

!local scalars
INTEGER :: clock, clock_rate, clock_start, clock_stop, i, j
CHARACTER(LEN=32) :: f
!local arrays
INTEGER, DIMENSION(4) :: iseed
!intrinsic procedures
INTRINSIC :: COUNT, DBLE, MAX, MAXVAL, MOD, NEW_LINE, SYSTEM_CLOCK
!external procedures
REAL(dp) :: dlange, dznrm2
EXTERNAL :: dlange, dznrm2
!file location (where problem files are stored)
CHARACTER(*), PARAMETER :: fileplace1="/home/thomas/Documents/FORTRAN/LMPEP2/PROBLEMS/REAL/"
!file location (where results are stored)
CHARACTER(*), PARAMETER :: fileplace2="/home/thomas/Documents/FORTRAN/QUADEIG/quadratic-eigensolver/RESULTS/"

!create iseed, used in dlarnv and dlagge
CALL SYSTEM_CLOCK(COUNT=clock)
CALL srand(clock)
DO i=1,4
  iseed(i)=MOD(irand(),4095)
ENDDO
IF(MOD(iseed(4),2)==0) THEN
  iseed(4)=iseed(4)+1
ENDIF

!random problem or stored problem?
PRINT*, 'input 1 for random problem, 2 for stored problem'
READ*, opt

!create or open problem
IF(opt==1) THEN
  !read in size
  PRINT*, 'input size'
  READ*, n
  ALLOCATE(work(2*n), x(n), A(n,n), B(n,n), C(n,n))
  CALL dlarnv(2, iseed, n, x)
  CALL dlagge(n, n, n-1, n-1, x, C, n, iseed, work, info)
  CALL dlarnv(2, iseed, n, x)
  CALL dlagge(n, n, n-1, n-1, x, B, n, iseed, work, info)
  CALL dlarnv(2, iseed, n, x)
  CALL dlagge(n, n, n-1, n-1, x, A, n, iseed, work, info)
  DEALLOCATE(work, x)
  d=2
ELSEIF(opt==2) THEN
  !read in file name
  PRINT*, 'input file name'
  READ*, f
  OPEN(UNIT=1,FILE=fileplace1//f)
  !read in size and degree from file
  READ(1,*) n
  READ(1,*) d
  !read in scalar polynomial
  ALLOCATE(A(n,n), B(n,n), C(n,n))
  READ(1,*) C
  READ(1,*) B
  READ(1,*) A
  CLOSE(UNIT=1)
ENDIF

ALLOCATE(alphar(2*n), alphai(2*n), beta(2*n))
ALLOCATE(VL(n,2*n),VR(n,2*n),s(2*n),beVL(2*n),beVR(2*n))
tol=n*eps
CALL SYSTEM_CLOCK(count_rate=clock_rate)
CALL SYSTEM_CLOCK(COUNT=clock_start)
CALL DG3EVX(1, 'V', 'V', 0, tol, n,                          &
                   A, n, B, n, C, n,                                     &
                   alphar, alphai, beta, VL, n, VR, n, s, beVL, beVR,    &
                   iwarn, info)
CALL SYSTEM_CLOCK(COUNT=clock_stop)

!compute backward error and condition number using posterrcond
ALLOCATE(er(2*n),ei(2*n),xr(n,2*n),xi(n,2*n),yr(n,2*n),yi(n,2*n))
i=1
DO WHILE(i<=2*n)
  IF(DABS(alphai(i))<eps) THEN
    er(i)=alphar(i)/beta(i); ei(i)=alphai(i)/beta(i)
    xr(:,i)=VR(:,i); xi(:,i)=0
    yr(:,i)=VL(:,i); yi(:,i)=0
    tol=dznrm2(n,DCMPLX(xr(:,i)),1)
    xr(:,i)=xr(:,i)/tol
    tol=dznrm2(n,DCMPLX(yr(:,i)),1)
    yr(:,i)=yr(:,i)/tol
    i=i+1
  ELSE 
    er(i)=alphar(i)/beta(i); ei(i)=alphai(i)/beta(i)
    er(i+1)=alphar(i)/beta(i); ei(i+1)=-alphai(i)/beta(i)
    xr(:,i)=VR(:,i); xi(:,i)=VR(:,i+1)
    xr(:,i+1)=VR(:,i); xi(:,i+1)=-VR(:,i+1)
    yr(:,i)=VL(:,i); yi(:,i)=VL(:,i+1)
    yr(:,i+1)=VL(:,i); yi(:,i+1)=-VL(:,i+1)
    tol=dznrm2(n,DCMPLX(xr(:,i),xi(:,i)),1)
    xr(:,i)=xr(:,i)/tol; xi(:,i)=xi(:,i)/tol
    xr(:,i+1)=xr(:,i+1)/tol; xi(:,i+1)=xi(:,i+1)/tol
    tol=dznrm2(n,DCMPLX(yr(:,i),yi(:,i)),1)
    yr(:,i)=yr(:,i)/tol; yi(:,i)=yi(:,i)/tol
    yr(:,i+1)=yr(:,i+1)/tol; yi(:,i+1)=yi(:,i+1)/tol
    i=i+2
  ENDIF
ENDDO

ALLOCATE(p(n,3*n),berr(2*n),cond(2*n),ferr(2*n),ncoeff(3))
p(1:n,1:n)=C
p(1:n,n+1:2*n)=B
p(1:n,2*n+1:3*n)=A
DO i=1,3
  ncoeff(i)=dlange('F',n,n,p(1,n*(i-1)+1),n,x)
ENDDO
CALL dposterrcond(p, xr, xi, yr, yi, er, ei, ncoeff, berr, cond, ferr, d, n)

!save results
WRITE(*,'(20G15.4)') 'QUADEIG TIME =', DBLE(clock_stop-clock_start)/DBLE(clock_rate)
WRITE(*,'(20G15.4)') 'MAX BERR      =', MAXVAL(berr)
WRITE(*,'(20G15.4)') 'MAX FERR      =', MAXVAL(ferr)
OPEN(UNIT=1,FILE=fileplace2//"EIGVALUES")
DO i=1,n*d
  WRITE(1,*) er(i), ei(i)
ENDDO
 CLOSE(UNIT=1)
OPEN(UNIT=1,FILE=fileplace2//"BERR-COND-FERR")
DO i=1,n*d
  WRITE(1,*) berr(i), cond(i), ferr(i)
ENDDO
 CLOSE(UNIT=1)
OPEN(UNIT=1,FILE=fileplace2//"R-EIGVECTORS")
DO j=1,n*d
  DO i=1,n
    WRITE(1,*) xr(i,j), xi(i,j)
  ENDDO
  WRITE(1,*) NEW_LINE('c')
ENDDO
 CLOSE(UNIT=1)
OPEN(UNIT=1,FILE=fileplace2//"L-EIGVECTORS")
DO j=1,n*d
  DO i=1,n
    WRITE(1,*) yr(i,j), yi(i,j)
  ENDDO
  WRITE(1,*) NEW_LINE('c')
ENDDO
 CLOSE(UNIT=1)

!deallocate
!DEALLOCATE(p,berr,cond,ferr,ncoeff)
!DEALLOCATE(er,ei,xr,xi,yr,yi)
!DEALLOCATE(alphar, alphai, beta)
!DEALLOCATE(VL,VR,s,beVL,beVR)

CONTAINS

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
  a=d*(d-1)*p(:,1:n)
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

IF(der==0) THEN
  a=p(1)
  DO k=2,d+1
    a=t*a+p(k)
  ENDDO
ELSEIF(der==1) THEN
  a=d*p(1)
  DO k=2,d
    a=t*a+(d-k+1)*p(k)
  ENDDO
ELSE
  a=d*(d-1)*p(1)
  DO k=2,d-1
    a=t*a+(d-k+1)*(d-k)*p(k)
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
!			SUBROUTINE DCMOD				*
!************************************************************************
! Compute the module of the complex number a+bi, while avoiding harmul  *
! overflow and underflow.						*
!************************************************************************
FUNCTION dcmod(a, b) RESULT(r)
IMPLICIT NONE
!scalar arguments
REAL(dp), INTENT(IN) :: a, b
!local scalars
REAL(dp) :: r

IF(DABS(a)<eps .AND. DABS(b)<eps) THEN
  r=zero
  RETURN
ENDIF

IF(DABS(a)<DABS(b)) THEN
  r=DABS(b)*DSQRT(1+(a/b)**2)
ELSE
  r=DABS(a)*DSQRT(1+(b/a)**2)
ENDIF
RETURN
END FUNCTION dcmod

END PROGRAM dgelmpep_driver
