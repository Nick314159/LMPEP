PROGRAM quadeig_driver
USE environment
USE dgeeam_subroutines
IMPLICIT NONE

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
  OPEN(UNIT=1,FILE=problemsDir//f)
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
OPEN(UNIT=1,FILE=resultsDir//"EIGVALUES")
DO i=1,n*d
  WRITE(1,*) er(i), ei(i)
ENDDO
 CLOSE(UNIT=1)
OPEN(UNIT=1,FILE=resultsDir//"BERR-COND-FERR")
DO i=1,n*d
  WRITE(1,*) berr(i), cond(i), ferr(i)
ENDDO
 CLOSE(UNIT=1)
OPEN(UNIT=1,FILE=resultsDir//"R-EIGVECTORS")
DO j=1,n*d
  DO i=1,n
    WRITE(1,*) xr(i,j), xi(i,j)
  ENDDO
  WRITE(1,*) NEW_LINE('c')
ENDDO
 CLOSE(UNIT=1)
OPEN(UNIT=1,FILE=resultsDir//"L-EIGVECTORS")
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
END PROGRAM quadeig_driver

