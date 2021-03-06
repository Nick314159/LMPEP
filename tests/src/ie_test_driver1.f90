PROGRAM ie_test_driver1
USE environment
USE dgelmpep_subroutines
IMPLICIT NONE

INTEGER :: maxSize, maxDegree, startingSize, startingDegree
CHARACTER(len=32) :: arg
!local scalars
INTEGER :: clock, clock_rate, clock_start, clock_stop, d, i, info, n, m, j
!local arrays
INTEGER, DIMENSION(4) :: iseed
REAL(dp), DIMENSION(:), ALLOCATABLE :: berr, cond, ferr, er, ei, ncoeff, work, x
REAL(dp), DIMENSION(:,:), ALLOCATABLE :: p, xr, xi, yr, yi
REAL(dp), DIMENSION(:,:), ALLOCATABLE :: timeStats
!intrinsic procedures
INTRINSIC :: COUNT, DBLE, MAX, MAXVAL, MOD, NEW_LINE, SYSTEM_CLOCK
!external procedures
REAL(dp) :: dlange
EXTERNAL :: dlange
!create iseed, used in dlarnv and dlagge
CALL SYSTEM_CLOCK(COUNT=clock)
CALL srand(clock)
DO i=1,4
  iseed(i)=MOD(irand(),4095)
ENDDO
IF(MOD(iseed(4),2)==0) THEN
  iseed(4)=iseed(4)+1
ENDIF

CALL GETARG(1, arg)
READ (arg,'(I10)') startingSize
CALL GETARG(2, arg)
READ (arg,'(I10)') maxSize
CALL GETARG(3, arg)
READ (arg,'(I10)') startingDegree
CALL GETARG(4, arg)
READ (arg,'(I10)') maxDegree
CALL GETARG(5, arg)
READ (arg,'(I10)') m

ALLOCATE(timeStats(m,2))

OPEN(UNIT=1,FILE=resultsDir//"outputIepoly1Size.csv")
WRITE(1, '(A)',  advance='no') 'DEGREE,     '
WRITE(1, '(A)',  advance='no') 'SIZE,    '
WRITE(1, '(A)',  advance='no') 'NR TIME,         '
WRITE(1, '(A)',  advance='no') 'NP TIME         '
WRITE(1, *)

  d = 2
  n = startingSize
  DO WHILE (n<maxSize)
    WRITE(1, '(i10)', advance='no') d
    WRITE(1, '(A)', advance='no') ', '
    WRITE(1, '(i10)', advance='no') n
    WRITE(1, '(A)', advance='no') ', '
    DO j=1,m
      !create problem
      ALLOCATE(work(n+n*(d+1)), x(n), p(n,n*(d+1)))
      DO i=1,d+1
        CALL dlarnv(2, iseed, n, x)
        CALL dlagge(n, n, n-1, n-1, x, p(1,n*(i-1)+1), n, iseed, work, info)
      ENDDO
      DEALLOCATE(work, x)
    
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
      !CALL dposterrcond(p, xr, xi, yr, yi, er, ei, ncoeff, berr, cond, ferr, d, n)
      !record results
      timeStats(j,1) = DBLE(clock_stop-clock_start)/DBLE(clock_rate)
    
      DEALLOCATE(xr, xi, yr, yi, berr, er, ei, ncoeff, cond, ferr)
    
     !solve problem using Laguerre's Method
      ALLOCATE(xr(n,n*d), xi(n,n*d), yr(n,n*d), yi(n,n*d))
      ALLOCATE(berr(n*d), er(n*d), ei(n*d), ncoeff(n*d), cond(n*d), ferr(n*d))
      DO i=1,d+1
        ncoeff(i)=dlange('F',n,n,p(1,n*(i-1)+1),n,x)
      ENDDO
      CALL SYSTEM_CLOCK(count_rate=clock_rate)
      CALL SYSTEM_CLOCK(COUNT=clock_start)
      CALL dgelm(p, xr, xi, yr, yi, er, ei, berr, ncoeff, iseed, d, n, 'NP')
      CALL SYSTEM_CLOCK(COUNT=clock_stop)
    
      !bacward error, condition number for Laguerre's Method
      !CALL dposterrcond(p, xr, xi, yr, yi, er, ei, ncoeff, berr, cond, ferr, d, n)
      !record results
      timeStats(j,2) = DBLE(clock_stop-clock_start)/DBLE(clock_rate)
    
      DEALLOCATE(p, xr, xi, yr, yi, berr, er, ei, ncoeff, cond, ferr)
    END DO
    !print results
    WRITE(1,'(20G15.4)', advance='no') SUM(timeStats(:,1))/m
    WRITE(1, '(A)', advance='no') ', '
    WRITE(1,'(20G15.4)', advance='no') SUM(timeStats(:,2))/m
    WRITE(1, *) 
    n = 2* n
  END DO
  CLOSE(UNIT=1) 
  
 
 
OPEN(UNIT=1,FILE=resultsDir//"outputIepoly1Degree.csv")
WRITE(1, '(A)',  advance='no') 'DEGREE,     '
WRITE(1, '(A)',  advance='no') 'SIZE,    '
WRITE(1, '(A)',  advance='no') 'NR TIME,         '
WRITE(1, '(A)',  advance='no') 'NP TIME         '
WRITE(1, *)

  d = startingDegree
  n = 2
  DO WHILE (d<maxDegree)
    WRITE(1, '(i6)', advance='no') d
    WRITE(1, '(A)', advance='no') ', '
    WRITE(1, '(i6)', advance='no') n
    WRITE(1, '(A)', advance='no') ', '
    DO j=1,m
      !create problem
      ALLOCATE(work(n+n*(d+1)), x(n), p(n,n*(d+1)))
      DO i=1,d+1
        CALL dlarnv(2, iseed, n, x)
        CALL dlagge(n, n, n-1, n-1, x, p(1,n*(i-1)+1), n, iseed, work, info)
      ENDDO
      DEALLOCATE(work, x)
    
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
      !CALL dposterrcond(p, xr, xi, yr, yi, er, ei, ncoeff, berr, cond, ferr, d, n)
      !record results
      timeStats(j,1) = DBLE(clock_stop-clock_start)/DBLE(clock_rate)
    
      DEALLOCATE(xr, xi, yr, yi, berr, er, ei, ncoeff, cond, ferr)
    
      !solve problem using Laguerre's Method
      ALLOCATE(xr(n,n*d), xi(n,n*d), yr(n,n*d), yi(n,n*d))
      ALLOCATE(berr(n*d), er(n*d), ei(n*d), ncoeff(n*d), cond(n*d), ferr(n*d))
      DO i=1,d+1
        ncoeff(i)=dlange('F',n,n,p(1,n*(i-1)+1),n,x)
      ENDDO
      CALL SYSTEM_CLOCK(count_rate=clock_rate)
      CALL SYSTEM_CLOCK(COUNT=clock_start)
      CALL dgelm(p, xr, xi, yr, yi, er, ei, berr, ncoeff, iseed, d, n, 'NP')
      CALL SYSTEM_CLOCK(COUNT=clock_stop)
    
      !bacward error, condition number for Laguerre's Method
      !CALL dposterrcond(p, xr, xi, yr, yi, er, ei, ncoeff, berr, cond, ferr, d, n)
      !record results
      timeStats(j,2) = DBLE(clock_stop-clock_start)/DBLE(clock_rate)
      
      DEALLOCATE(p, xr, xi, yr, yi, berr, er, ei, ncoeff, cond, ferr)
    END DO
    !print results
    WRITE(1,'(20G15.4)', advance='no') SUM(timeStats(:,1))/m
    WRITE(1, '(A)', advance='no') ', '
    WRITE(1,'(20G15.4)', advance='no') SUM(timeStats(:,2))/m
    WRITE(1, *) 
    d = 2 * d
  END DO
  CLOSE(UNIT=1) 
   
  
END PROGRAM ie_test_driver1
