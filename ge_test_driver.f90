PROGRAM ge_test_driver
USE dgeeam_subroutines
IMPLICIT NONE

!local scalars
INTEGER :: c, clock, clock_rate, clock_start, clock_stop, d, i, info, j, n
CHARACTER (LEN=32) f
!local arrays
INTEGER, DIMENSION(4) :: iseed
REAL(dp), DIMENSION(:), ALLOCATABLE :: berr, cond, ferr, er, ei, ncoeff, work, x
REAL(dp), DIMENSION(:,:), ALLOCATABLE :: p, xr, xi, yr, yi
!intrinsic procedures
INTRINSIC :: COUNT, DBLE, MAX, MAXVAL, MOD, NEW_LINE, SYSTEM_CLOCK
!external procedures
REAL(dp) :: dlange
EXTERNAL :: dlange
!file location (where problem files are stored)
CHARACTER(*), PARAMETER :: fileplace1="/home/thomas/Documents/FORTRAN/Nick/LMPEPtests/PROBLEMS/REAL/"
!file location (where results are stored)
CHARACTER(*), PARAMETER :: fileplace2="/home/thomas/Documents/FORTRAN/Nick/LMPEPtests/results/"

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
READ*, c

!create or open problem
IF(c==1) THEN
  !read in degree
  PRINT*, 'input size and degree'
  READ*, n
  READ*, d
  ALLOCATE(work(n+n*(d+1)), x(n), p(n,n*(d+1)))
  DO i=1,d+1
    CALL dlarnv(2, iseed, n, x)
    CALL dlagge(n, n, n-1, n-1, x, p(1,n*(i-1)+1), n, iseed, work, info)
  ENDDO
  DEALLOCATE(work, x)
ELSEIF(c==2) THEN
  !read in file name
  PRINT*, 'input file name'
  READ*, f
  OPEN(UNIT=1,FILE=fileplace1//f)
  !read in size and degree from file
  READ(1,*) n
  READ(1,*) d
  !read in scalar polynomial
  ALLOCATE(p(n,n*(d+1)))
  READ(1,*) p
  CLOSE(UNIT=1)
ENDIF

!solve problem using Laguerre's Method
ALLOCATE(xr(n,n*d), xi(n,n*d), yr(n,n*d), yi(n,n*d))
ALLOCATE(berr(n*d), er(n*d), ei(n*d), ncoeff(n*d), cond(n*d), ferr(n*d))
DO i=1,d+1
  ncoeff(i)=dlange('F',n,n,p(1,n*(i-1)+1),n,x)
ENDDO
CALL SYSTEM_CLOCK(count_rate=clock_rate)
CALL SYSTEM_CLOCK(COUNT=clock_start)
CALL dgelm(p, xr, xi, yr, yi, er, ei, berr, ncoeff, iseed, d, n, 'NR')
CALL SYSTEM_CLOCK(COUNT=clock_stop)

!bacward error, condition number for Laguerre's Method
CALL dposterrcond(p, xr, xi, yr, yi, er, ei, ncoeff, berr, cond, ferr, d, n)
!print results
WRITE(*,'(20G15.4)') 'Laguerre TIME =', DBLE(clock_stop-clock_start)/DBLE(clock_rate)
WRITE(*,'(20G15.4)') 'MAX BERR      =', MAXVAL(berr)
WRITE(*,'(20G15.4)') 'MAX FERR      =', MAXVAL(ferr)
!newline
WRITE(*,*) NEW_LINE('c')
!solve problem using Ehrlich-Aberth method
CALL SYSTEM_CLOCK(count_rate=clock_rate)
CALL SYSTEM_CLOCK(COUNT=clock_start)
CALL dgeeam(p, xr, xi, yr, yi, er, ei, berr, ncoeff, iseed, d, n, 'NP')
CALL SYSTEM_CLOCK(COUNT=clock_stop)

!bacward error, condition number for Ehrlich-Aberth method
CALL dposterrcond(p, xr, xi, yr, yi, er, ei, ncoeff, berr, cond, ferr, d, n)
!print results
WRITE(*,'(20G15.4)') 'Ehrlich-Aberth TIME =', DBLE(clock_stop-clock_start)/DBLE(clock_rate)
WRITE(*,'(20G15.4)') 'MAX BERR            =', MAXVAL(berr)
WRITE(*,'(20G15.4)') 'MAX FERR            =', MAXVAL(ferr)







END PROGRAM ge_test_driver
