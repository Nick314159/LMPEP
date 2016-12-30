PROGRAM ie_test_driver2
USE dgeeam_subroutines
IMPLICIT NONE

!local scalars
INTEGER :: c, clock, clock_rate, clock_start, clock_stop, d, i, info, j, n, k
CHARACTER (LEN=32) f
!local arrays
INTEGER, DIMENSION(4) :: iseed
REAL(dp), DIMENSION(:), ALLOCATABLE :: berr, cond, ferr, er, ei, ncoeff, work, x, ier, iei
REAL(dp), DIMENSION(:,:), ALLOCATABLE :: p, xr, xi, yr, yi
!intrinsic procedures
INTRINSIC :: COUNT, DBLE, MAX, MAXVAL, MOD, NEW_LINE, SYSTEM_CLOCK
!external procedures
REAL(dp) :: dlange
EXTERNAL :: dlange
!file location (where problem files are stored)
CHARACTER(*), PARAMETER :: fileplace1="/home/thomas/Documents/FORTRAN/Nick/LMPEPtests/PROBLEMS/REAL/"
!CHARACTER(*), PARAMETER :: fileplace1="/home/nsteckley/Documents/Personal/Cameron/LMPEP/tests/PROBLEMS/REAL/"
!file location (where results are stored)
CHARACTER(*), PARAMETER :: fileplace2="/home/thomas/Documents/FORTRAN/Nick/LMPEPtests/results/"
!CHARACTER(*), PARAMETER :: fileplace2="/home/nsteckley/Documents/Personal/Cameron/LMPEP/tests/results/"
CHARACTER (LEN=64), DIMENSION(4) :: tests

!create iseed, used in dlarnv and dlagge
CALL SYSTEM_CLOCK(COUNT=clock)
CALL srand(clock)
DO i=1,4
  iseed(i)=MOD(irand(),4095)
ENDDO
IF(MOD(iseed(4),2)==0) THEN
  iseed(4)=iseed(4)+1
ENDIF
tests(1) = 'bilby.txt'
tests(2) = 'butterfly.txt'
tests(3) = 'cd_player.txt'
tests(4) = 'spring.txt'

OPEN(UNIT=2, FILE=fileplace2//'outputIepoly2.txt')
DO k=1,4
  OPEN(UNIT=1,FILE=fileplace1//tests(k))
  !read in size and degree from file
  READ(1,*) n
  READ(1,*) d
  !read in scalar polynomial
  ALLOCATE(p(n,n*(d+1)))
  READ(1,*) p
  CLOSE(UNIT=1)

  !solve problem using Laguerre's Method
  ALLOCATE(xr(n,n*d), xi(n,n*d), yr(n,n*d), yi(n,n*d))
  ALLOCATE(berr(n*d), er(n*d), ei(n*d), ncoeff(n*d), cond(n*d), ferr(n*d),  ier(n*d), iei(n*d))
  DO i=1,d+1
    ncoeff(i)=dlange('F',n,n,p(1,n*(i-1)+1),n,x)
  ENDDO
  CALL SYSTEM_CLOCK(count_rate=clock_rate)
  CALL SYSTEM_CLOCK(COUNT=clock_start)
  CALL dgelmt(p, xr, xi, yr, yi, er, ei, berr, ncoeff, iseed, d, n, 'NR', ier, iei)
  CALL SYSTEM_CLOCK(COUNT=clock_stop)

  !bacward error, condition number for Laguerre's Method
  CALL dposterrcond(p, xr, xi, yr, yi, er, ei, ncoeff, berr, cond, ferr, d, n)
  !print results
  
  WRITE(2,*) tests(k)
  WRITE(2,*) 'NR'
  WRITE(2,*) ier
  WRITE(2,*)
  WRITE(2,*) iei
  WRITE(2,*)
  WRITE(2,*) er
  WRITE(2,*)
  WRITE(2,*) ei
  WRITE(2,*)

  DEALLOCATE(xr, xi, yr, yi, berr, er, ei, ncoeff, cond, ferr, ier, iei)
  
   !solve problem using Laguerre's Method
  ALLOCATE(xr(n,n*d), xi(n,n*d), yr(n,n*d), yi(n,n*d))
  ALLOCATE(berr(n*d), er(n*d), ei(n*d), ncoeff(n*d), cond(n*d), ferr(n*d),  ier(n*d), iei(n*d))
  DO i=1,d+1
    ncoeff(i)=dlange('F',n,n,p(1,n*(i-1)+1),n,x)
  ENDDO
  CALL SYSTEM_CLOCK(count_rate=clock_rate)
  CALL SYSTEM_CLOCK(COUNT=clock_start)
  CALL dgelmt(p, xr, xi, yr, yi, er, ei, berr, ncoeff, iseed, d, n, 'NP', ier, iei)
  CALL SYSTEM_CLOCK(COUNT=clock_stop)

  !bacward error, condition number for Laguerre's Method
  CALL dposterrcond(p, xr, xi, yr, yi, er, ei, ncoeff, berr, cond, ferr, d, n)
  !print results
  
  WRITE(2,*) tests(k)
  WRITE(2,*) 'NP'
  WRITE(2,*) ier
  WRITE(2,*)
  WRITE(2,*) iei
  WRITE(2,*)
  WRITE(2,*) er
  WRITE(2,*)
  WRITE(2,*) ei
  WRITE(2,*)

  DEALLOCATE(p, xr, xi, yr, yi, berr, er, ei, ncoeff, cond, ferr, ier, iei)
END DO 
 CLOSE(UNIT=2)
END PROGRAM ie_test_driver2
