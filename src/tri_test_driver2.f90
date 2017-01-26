PROGRAM tri_test_driver2
USE environment
USE dgtlmpep_subroutines
USE eigensolve
IMPLICIT NONE

!===Variables===
!LMPEP
INTEGER :: clock, clock_rate, clock_start, clock_stop
INTEGER, DIMENSION(4) :: iseed
REAL(dp), DIMENSION(:), ALLOCATABLE :: er, ei, ncoeff, berr, ferr
REAL(dp), DIMENSION(:,:), ALLOCATABLE :: pd, pdl, pdu, xr, xi, yr, yi
!EIGENSOLVE
INTEGER                  :: n,i,j,m, jmax, jmin
REAL(dp)                 :: alpha
REAL(dp),DIMENSION(:)    :: a,s,rad ,ss,tt,cond
COMPLEX(dp),DIMENSION(:) :: z, co, si, ad, adl, adu, x, y
ALLOCATABLE              :: a, s, z, rad,ss,tt,cond, co, si, ad, adl, adu, x, y
REAL(dp)                 :: iter ,aaa,bbb,theta,h
CHARACTER(len=20)        :: filename
REAL                     :: ru
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

!!! Accuracy Tests
!size
n=100
DO j=1,9
  PRINT*, 'test ', j
!!! Allocate
  ALLOCATE(a(n), s(n), z(n), cond(n), rad(n))
  ALLOCATE(pdl(n-1,2), pd(n,2), pdu(n-1,2))
  ALLOCATE(ncoeff(2))
  ALLOCATE(xr(n,n), xi(n,n), yr(n,n), yi(n,n))
  ALLOCATE(berr(n), ferr(n), er(n), ei(n))
  IF(j==1)THEN
!!! TEST 1
    DO i=1,n
      a(i)=i*(-1.d0)**(i/8)
      s(i)=((-1.d0)**i)/i
    END DO
  END IF
  IF(j==2)THEN
!!! TEST 2
    DO i=1,n
      a(i)= 10*(-1.d0)**(i/8)
      s(i)=i*(-1)**(i/9)
    END DO
  END IF
  IF(j==3)THEN
!!! TEST 3
    DO i=1,n
      a(i)=i
      s(i)=(n-i+1)
    END DO
  END IF
  IF(j==4)THEN
!!! TEST 4
    DO i=1,n
      a(i)=1.d0*(-1)**i
      s(i)=20.d0*(-1)**(i/5)
    END DO
  END IF
  IF(j==5)THEN
!!! TEST 5. 
    DO i=1,n
      a(i)=(-1)**(i/4)*10.d0**(5*(-1)**i)
      s(i)=(-1)**(i/3)
    END DO
  END IF
  IF(j==6)THEN
!!! TEST 6
    DO i=1,n
      a(i)=2
      s(i)=1
    END DO
  END IF
  IF(j==7)THEN
!!! TEST 7
    DO i=1,n
      a(i)=1.d0/i+1.d0/(n-i+1)
      s(i)=(1.d0/i)*(-1)**(i/9)
    END DO
  END IF
  IF(j==8)THEN
!!! TEST 8
    DO i=1,n
      a(i)=i*(-1)**(i/13)*(-1)**(i/5)
      s(i)=(-1)**(i/11)*(n-i+1)**2.d0
    END DO
  END IF
  IF(j==9)THEN
!!! * TEST 9 The matrix {{-x,1,0},{1,x,1},{0,1,x}} is such that det =-x^3
!!!  the matrix obtained by concatenating two blocks equal to the previous
!!!  matrix has nonzero eigenvalues 
    DO i=1,n
      a(i)=1
      s(i)=1.d0
      IF(i>=n/2)s(i)=-1.d0
    END DO
  END IF
!!! Compute eigenvalues using eigen
  CALL eigen(n,a,s,z,cond)
!!! Error estimates for eigen
  pdl(:,1)=one; pd(:,1)=a; pdu(:,1)=one
  pdl(:,2)=zero; pd(:,2)=-s; pdu(:,2)=zero
  ALLOCATE(ad(n), adu(n-1), adl(n-1), co(n-1), si(n-1), x(n), y(n))
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
  DEALLOCATE(adl, ad, adu, co, si, x, y)
  !post err analysis
  CALL dposterrcond(pdl, pd, pdu, xr, xi, yr, yi, er, ei, ncoeff, berr, cond, ferr, 1, n)
  PRINT*, 'MAX FERR =', MAXVAL(ferr)
!!! Deallocate
  DEALLOCATE(berr, ferr, er, ei)
  DEALLOCATE(xr, xi, yr, yi)
  DEALLOCATE(ncoeff)
  DEALLOCATE(pdl, pd, pdu)
  DEALLOCATE(z, cond, rad)
!!! Allocate
  ALLOCATE(z(n), cond(n), rad(n))
  ALLOCATE(pdl(n-1,2), pd(n,2), pdu(n-1,2))
  ALLOCATE(ncoeff(2))
  ALLOCATE(xr(n,n), xi(n,n), yr(n,n), yi(n,n))
  ALLOCATE(berr(n), ferr(n), er(n), ei(n))
!!! Compute eigenvalues using DGTLMPEP
  !coefficient norm
  DO i=1,2
    ncoeff(i)=dlangt('F',n,pdl(1,i),pd(1,i),pdu(1,i))
  ENDDO
  !solve problem
  CALL dgtlm(pdl, pd, pdu, xr, xi, yr, yi, er, ei, berr, ncoeff, iseed, 1, n)
  CALL dposterrcond(pdl, pd, pdu, xr, xi, yr, yi, er, ei, ncoeff, berr, cond, ferr, 1, n)
  PRINT*, 'MAX FERR =', MAXVAL(ferr)
!!! Deallocate
  DEALLOCATE(berr, ferr, er, ei)
  DEALLOCATE(xr, xi, yr, yi)
  DEALLOCATE(ncoeff)
  DEALLOCATE(pdl, pd, pdu)
  DEALLOCATE(a, s, z, cond, rad)
ENDDO


END PROGRAM tri_test_driver2
