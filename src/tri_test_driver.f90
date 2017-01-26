PROGRAM tri_test_driver
USE dgtlmpep_subroutines
USE eigensolve
IMPLICIT NONE

!===Variables===
!LMPEP
INTEGER :: clock, clock_rate, clock_start, clock_stop, d
INTEGER, DIMENSION(4) :: iseed
REAL(dp), DIMENSION(:), ALLOCATABLE :: er, ei, ncoeff, berr, ferr
REAL(dp), DIMENSION(:,:), ALLOCATABLE :: pd, pdl, pdu, xr, xi, yr, yi
!EIGENSOLVE
INTEGER                              :: n,i,j, jmax, jmin
REAL(dp),DIMENSION(:),ALLOCATABLE    :: a,b,c,s,cond
COMPLEX(dp),DIMENSION(:),ALLOCATABLE :: z, co, si, ad, adl, adu, x, y
REAL(dp)                             :: alpha
!external procedures
REAL(dp) :: dlangt
EXTERNAL :: dlangt

!create iseed, used in dgtlmpep
CALL SYSTEM_CLOCK(COUNT=clock)
CALL srand(clock)
DO i=1,4
  iseed(i)=MOD(irand(),4095)
ENDDO
IF(MOD(iseed(4),2)==0) THEN
  iseed(4)=iseed(4)+1
ENDIF

!===DGTLMPEP Accuracy Tests===
PRINT*, 'DGTLMPEP Accuracy Tests'
n=20
d=1
DO j=1,9
  PRINT*, 'test ', j
!!! Allocate
  ALLOCATE(a(n), s(n))
  ALLOCATE(pdl(n-1,d+1), pd(n,d+1), pdu(n-1,d+1))
  ALLOCATE(ncoeff(d+1))
  ALLOCATE(xr(n,n*d), xi(n,n*d), yr(n,n*d), yi(n,n*d))
  ALLOCATE(berr(n*d), cond(n*d), ferr(n*d), er(n*d), ei(n*d))
  IF(j==1)THEN
!!! TEST 1
    DO i=1,n
      a(i)=i*(-1.d0)**(i/8)
      s(i)=((-1.d0)**i)/i
    ENDDO
  ENDIF
  IF(j==2)THEN
!!! TEST 2
    DO i=1,n
      a(i)= 10*(-1.d0)**(i/8)
      s(i)=i*(-1)**(i/9)
    ENDDO
  ENDIF
  IF(j==3)THEN
!!! TEST 3
    DO i=1,n
      a(i)=i
      s(i)=(n-i+1)
    ENDDO
  ENDIF
  IF(j==4)THEN
!!! TEST 4
    DO i=1,n
      a(i)=1.d0*(-1)**i
      s(i)=20.d0*(-1)**(i/5)
    ENDDO
  ENDIF
  IF(j==5)THEN
!!! TEST 5. 
    DO i=1,n
      a(i)=(-1)**(i/4)*10.d0**(5*(-1)**i)
      s(i)=(-1)**(i/3)
    ENDDO
  ENDIF
  IF(j==6)THEN
!!! TEST 6
    DO i=1,n
      a(i)=2
      s(i)=1
    ENDDO
  ENDIF
  IF(j==7)THEN
!!! TEST 7
    DO i=1,n
      a(i)=1.d0/i+1.d0/(n-i+1)
      s(i)=(1.d0/i)*(-1)**(i/9)
    ENDDO
  ENDIF
  IF(j==8)THEN
!!! TEST 8
    DO i=1,n
      a(i)=i*(-1)**(i/13)*(-1)**(i/5)
      s(i)=(-1)**(i/11)*(n-i+1)**2.d0
    ENDDO
  ENDIF
  IF(j==9)THEN
!!! * TEST 9 The matrix {{-x,1,0},{1,x,1},{0,1,x}} is such that det =-x^3
!!!  the matrix obtained by concatenating two blocks equal to the previous
!!!  matrix has nonzero eigenvalues 
    DO i=1,n
      a(i)=1
      s(i)=1.d0
      IF(i>=n/2)s(i)=-1.d0
    ENDDO
  ENDIF

!!! Compute eigenvalues using DGTLMPEP
  pdl(:,1)=one; pd(:,1)=a; pdu(:,1)=one
  pdl(:,2)=zero; pd(:,2)=-s; pdu(:,2)=zero
  !coefficient norm
  DO i=1,d+1
    ncoeff(i)=dlangt('F',n,pdl(1,i),pd(1,i),pdu(1,i))
  ENDDO
  !solve problem
  CALL dgtlm(pdl,pd,pdu,xr,xi,yr,yi,er,ei,berr,ncoeff,iseed,d,n)
  !error estimates for dgtlmpep
  CALL dposterrcond(pdl,pd,pdu,xr,xi,yr,yi,er,ei,ncoeff,berr,cond,ferr,d,n)
  PRINT*, 'MAX FERR =', MAXVAL(ferr)

!!! Deallocate
  DEALLOCATE(a, s)
  DEALLOCATE(pdl, pd, pdu)
  DEALLOCATE(ncoeff)
  DEALLOCATE(xr, xi, yr, yi)
  DEALLOCATE(berr, cond, ferr, er, ei)
ENDDO

PRINT*, 

!===EIGEN Accuracy Tests===
PRINT*, 'EIGEN Accuracy Tests'
n=100
d=1
DO j=1,9
  PRINT*, 'test ', j
!!! Allocate
  ALLOCATE(a(n), s(n), z(n))
  ALLOCATE(pdl(n-1,d+1), pd(n,d+1), pdu(n-1,d+1))
  ALLOCATE(ncoeff(d+1))
  ALLOCATE(xr(n,n*d), xi(n,n*d), yr(n,n*d), yi(n,n*d))
  ALLOCATE(berr(n*d), cond(n*d), ferr(n*d), er(n*d), ei(n*d))
  ALLOCATE(ad(n), adu(n-1), adl(n-1), co(n-1), si(n-1), x(n), y(n))
  IF(j==1)THEN
!!! TEST 1
    DO i=1,n
      a(i)=i*(-1.d0)**(i/8)
      s(i)=((-1.d0)**i)/i
    ENDDO
  ENDIF
  IF(j==2)THEN
!!! TEST 2
    DO i=1,n
      a(i)= 10*(-1.d0)**(i/8)
      s(i)=i*(-1)**(i/9)
    ENDDO
  ENDIF
  IF(j==3)THEN
!!! TEST 3
    DO i=1,n
      a(i)=i
      s(i)=(n-i+1)
    ENDDO
  ENDIF
  IF(j==4)THEN
!!! TEST 4
    DO i=1,n
      a(i)=1.d0*(-1)**i
      s(i)=20.d0*(-1)**(i/5)
    ENDDO
  ENDIF
  IF(j==5)THEN
!!! TEST 5. 
    DO i=1,n
      a(i)=(-1)**(i/4)*10.d0**(5*(-1)**i)
      s(i)=(-1)**(i/3)
    ENDDO
  ENDIF
  IF(j==6)THEN
!!! TEST 6
    DO i=1,n
      a(i)=2
      s(i)=1
    ENDDO
  ENDIF
  IF(j==7)THEN
!!! TEST 7
    DO i=1,n
      a(i)=1.d0/i+1.d0/(n-i+1)
      s(i)=(1.d0/i)*(-1)**(i/9)
    ENDDO
  ENDIF
  IF(j==8)THEN
!!! TEST 8
    DO i=1,n
      a(i)=i*(-1)**(i/13)*(-1)**(i/5)
      s(i)=(-1)**(i/11)*(n-i+1)**2.d0
    ENDDO
  ENDIF
  IF(j==9)THEN
!!! * TEST 9 The matrix {{-x,1,0},{1,x,1},{0,1,x}} is such that det =-x^3
!!!  the matrix obtained by concatenating two blocks equal to the previous
!!!  matrix has nonzero eigenvalues 
    DO i=1,n
      a(i)=1
      s(i)=1.d0
      IF(i>=n/2)s(i)=-1.d0
    ENDDO
  ENDIF

!!! Compute eigenvalues using EIGEN
  CALL eigen(n,a,s,z,cond)
  !error estimates for eigen
  pdl(:,1)=one; pd(:,1)=a; pdu(:,1)=one
  pdl(:,2)=zero; pd(:,2)=-s; pdu(:,2)=zero
  !coefficient norm
  DO i=1,d+1
    ncoeff(i)=dlangt('F',n,pdl(1,i),pd(1,i),pdu(1,i))
  ENDDO
  !eigenvectors
  DO i=1,n
    er(i)=DBLE(z(i)); ei(i)=DIMAG(z(i))
    IF(ZABS(z(i))>one) THEN
      z(i)=1/z(i)
      CALL zrevgteval(pdl, pd, pdu, z(i), adl, ad, adu, 1, n, 0)
      CALL drevseval(ncoeff, ZABS(z(i)), alpha, 1, 0)
      adl=adl/alpha; ad=ad/alpha; adu=adu/alpha
      CALL zgtqr(adl, ad, adu, co, si, n)
      jmax=ZGTJMAX(ad,n)
      jmin=ZGTJMIN(ad,n)
      IF(ZABS(ad(jmin))<eps) THEN
        CALL zgtker1(adl, ad, adu, x, y, co, si, jmin, jmax, n)
        xr(1:n,i)=DBLE(x); xi(1:n,i)=DIMAG(x)
        yr(1:n,i)=DBLE(y); yi(1:n,i)=DIMAG(y)
      ELSE
        CALL zgtker2(adl, ad, adu, x, y, co, si, jmin, jmax, n)
        xr(1:n,i)=DBLE(x); xi(1:n,i)=DIMAG(x)
        yr(1:n,i)=DBLE(y); yi(1:n,i)=DIMAG(y)
      ENDIF
    ELSE
      CALL zgteval(pdl, pd, pdu, z(i), adl, ad, adu, 1, n, 0)
      CALL dseval(ncoeff, ZABS(z(i)), alpha, 1, 0)
      adl=adl/alpha; ad=ad/alpha; adu=adu/alpha
      CALL zgtqr(adl, ad, adu, co, si, n)
      jmax=ZGTJMAX(ad,n)
      jmin=ZGTJMIN(ad,n)
      IF(ZABS(ad(jmin))<eps) THEN
        CALL zgtker1(adl, ad, adu, x, y, co, si, jmin, jmax, n)
        xr(1:n,i)=DBLE(x); xi(1:n,i)=DIMAG(x)
        yr(1:n,i)=DBLE(y); yi(1:n,i)=DIMAG(y)
      ELSE
        CALL zgtker2(adl, ad, adu, x, y, co, si, jmin, jmax, n)
        xr(1:n,i)=DBLE(x); xi(1:n,i)=DIMAG(x)
        yr(1:n,i)=DBLE(y); yi(1:n,i)=DIMAG(y)
      ENDIF
    ENDIF
  ENDDO
  !error estimates for dgtlmpep
  CALL dposterrcond(pdl,pd,pdu,xr,xi,yr,yi,er,ei,ncoeff,berr,cond,ferr,d,n)
  PRINT*, 'MAX FERR =', MAXVAL(ferr)

!!! Deallocate
  DEALLOCATE(a, s, z)
  DEALLOCATE(pdl, pd, pdu)
  DEALLOCATE(ncoeff)
  DEALLOCATE(xr, xi, yr, yi)
  DEALLOCATE(berr, cond, ferr, er, ei)
  DEALLOCATE(ad, adu, adl, co, si, x, y)
ENDDO

PRINT*,

!===Cost Complexity Tests===
PRINT*, 'Cost Complexity Tests'
n=10
d=1
DO WHILE(n<200)
PRINT*, 'size ', n
!!! Allocate
  ALLOCATE(a(n), b(n-1), c(n-1), s(n), z(n))
  ALLOCATE(pdl(n-1,d+1), pd(n,d+1), pdu(n-1,d+1))
  ALLOCATE(ncoeff(d+1))
  ALLOCATE(xr(n,n*d), xi(n,n*d), yr(n,n*d), yi(n,n*d))
  ALLOCATE(berr(n*d), cond(n*d), er(n*d), ei(n*d))
!!! Create random problem
  !DO i=1,n
  !  CALL RANDOM_NUMBER(ru)
  !  a(i)=0.5-ru
  !  CALL RANDOM_NUMBER(ru)
  !  s(i)=0.5-ru
  !ENDDO
  CALL dlarnv(2, iseed, n, a)
  CALL dlarnv(2, iseed, n-1, b)
  CALL dlarnv(2, iseed, n-1, c)
  CALL dlarnv(2, iseed, n, s)
!!! Compute eigenvalues using DGTLMPEP
  pdl(:,1)=c; pd(:,1)=a; pdu(:,1)=b
  pdl(:,2)=zero; pd(:,2)=-s; pdu(:,2)=zero
  !coefficient norm
  DO i=1,d+1
    ncoeff(i)=dlangt('F',n,pdl(1,i),pd(1,i),pdu(1,i))
  ENDDO
  !solve problem
  CALL SYSTEM_CLOCK(count_rate=clock_rate)
  CALL SYSTEM_CLOCK(COUNT=clock_start)
  CALL dgtlm(pdl,pd,pdu,xr,xi,yr,yi,er,ei,berr,ncoeff,iseed,d,n)
  CALL SYSTEM_CLOCK(COUNT=clock_stop)
  PRINT*, 'DGTLMPEP time =', DBLE(clock_stop-clock_start)/DBLE(clock_rate)

!!! Compute eigenvalues using EIGEN
  CALL SYSTEM_CLOCK(count_rate=clock_rate)
  CALL SYSTEM_CLOCK(COUNT=clock_start)
  CALL normalize(n,a,b,c,s)
  CALL eigen(n,a,s,z,cond)
  CALL SYSTEM_CLOCK(COUNT=clock_stop)
  PRINT*, 'EIGEN time =', DBLE(clock_stop-clock_start)/DBLE(clock_rate)

!!! Deallocate
  DEALLOCATE(a, b, c, s, z)
  DEALLOCATE(pdl, pd, pdu)
  DEALLOCATE(ncoeff)
  DEALLOCATE(xr, xi, yr, yi)
  DEALLOCATE(berr, cond, er, ei)

!update n
n=2*n
ENDDO

END PROGRAM tri_test_driver
