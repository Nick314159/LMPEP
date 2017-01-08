PROGRAM cc_test_driver
USE environment
USE dhslmpep_subroutines
IMPLICIT NONE

INTEGER :: maxSize, maxDegree, startingSize, startingDegree
CHARACTER(LEN=32) :: arg
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

m=3
ALLOCATE(timeStats(m,2))

OPEN(UNIT=1,FILE=resultsDir//"outputComplexitySize.csv")
WRITE(1, '(A)',  advance='no') 'DEGREE,     '
WRITE(1, '(A)',  advance='no') 'SIZE,    '
WRITE(1, '(A)',  advance='no') 'GE TIME,          '
WRITE(1, '(A)',  advance='no') 'HS TIME         '
WRITE(1, *)


d = 2
n = startingSize
DO WHILE (n<maxSize)
  WRITE(1, '(i6)', advance='no') d
  WRITE(1, '(A)', advance='no') ', '
  WRITE(1, '(i6)', advance='no') n
  WRITE(1, '(A)', advance='no') ', '
    
  DO j=1,m
    !create general problem
    ALLOCATE(work(n+n*(d+1)), x(n), p(n,n*(d+1)), ncoeff(d+1))
    DO i=1,d+1
      CALL dlarnv(2, iseed, n, x)
      CALL dlagge(n, n, n-1, n-1, x, p(1,n*(i-1)+1), n, iseed, work, info)
    ENDDO
    DO i=1,d+1
      ncoeff(i)=dlange('F',n,n,p(1,n*(i-1)+1),n,x)
    ENDDO
    
    !solve problem using General code
    ALLOCATE(xr(n,n*d), xi(n,n*d), yr(n,n*d), yi(n,n*d))
    ALLOCATE(berr(n*d), er(n*d), ei(n*d), cond(n*d), ferr(n*d))
    CALL SYSTEM_CLOCK(count_rate=clock_rate)
    CALL SYSTEM_CLOCK(COUNT=clock_start)
    CALL dgelm(p, xr, xi, yr, yi, er, ei, berr, ncoeff, iseed, d, n, 'NR')
    CALL SYSTEM_CLOCK(COUNT=clock_stop)
    
    timeStats(j,1) = DBLE(clock_stop-clock_start)/DBLE(clock_rate)

    !create Hessenberg problem
    DO i=1,d+1
      CALL dlarnv(2, iseed, n, x)
      CALL dlagge(n, n, 1, n-1, x, p(1,n*(i-1)+1), n, iseed, work, info)
    ENDDO
    DO i=1,d+1
      ncoeff(i)=dlange('F',n,n,p(1,n*(i-1)+1),n,x)
    ENDDO
    DEALLOCATE(work, x)

    !solve problem using Hessenberg code
    CALL SYSTEM_CLOCK(count_rate=clock_rate)
    CALL SYSTEM_CLOCK(COUNT=clock_start)
    CALL dhslm(p, xr, xi, yr, yi, er, ei, berr, ncoeff, iseed, d, n)
    CALL SYSTEM_CLOCK(COUNT=clock_stop)
    
    timeStats(j,2) = DBLE(clock_stop-clock_start)/DBLE(clock_rate)
    
    DEALLOCATE(p, xr, xi, yr, yi, berr, er, ei, ncoeff, cond, ferr)
 
  END DO
  n = 2* n
  !print results
  WRITE(1,'(20G15.4)', advance='no') SUM(timeStats(:,1))/m
  WRITE(1, '(A)', advance='no') ', '
  WRITE(1,'(20G15.4)', advance='no') SUM(timeStats(:,2))/m
  WRITE(1, *) 
END DO
 CLOSE(UNIT=1) 
  
OPEN(UNIT=1,FILE=resultsDir//"outputComplexityDegree.csv")
WRITE(1, '(A)',  advance='no') 'DEGREE,     '
WRITE(1, '(A)',  advance='no') 'SIZE,    '
WRITE(1, '(A)',  advance='no') 'GE TIME,          '
WRITE(1, '(A)',  advance='no') 'HS TIME,         '
WRITE(1, *) 
  
n = 2
d = startingDegree
DO WHILE (d<maxDegree)
  WRITE(1, '(i6)', advance='no') d
  WRITE(1, '(A)', advance='no') ', '
  WRITE(1, '(i6)', advance='no') n
  WRITE(1, '(A)', advance='no') ', '
    DO j=1,m
    !create general problem
    ALLOCATE(work(n+n*(d+1)), x(n), p(n,n*(d+1)), ncoeff(d+1))
    DO i=1,d+1
      CALL dlarnv(2, iseed, n, x)
      CALL dlagge(n, n, n-1, n-1, x, p(1,n*(i-1)+1), n, iseed, work, info)
    ENDDO
    DO i=1,d+1
      ncoeff(i)=dlange('F',n,n,p(1,n*(i-1)+1),n,x)
    ENDDO
    
    !solve problem using General code
    ALLOCATE(xr(n,n*d), xi(n,n*d), yr(n,n*d), yi(n,n*d))
    ALLOCATE(berr(n*d), er(n*d), ei(n*d), cond(n*d), ferr(n*d))
    CALL SYSTEM_CLOCK(count_rate=clock_rate)
    CALL SYSTEM_CLOCK(COUNT=clock_start)
    CALL dgelm(p, xr, xi, yr, yi, er, ei, berr, ncoeff, iseed, d, n, 'NR')
    CALL SYSTEM_CLOCK(COUNT=clock_stop)
    
    timeStats(j,1) = DBLE(clock_stop-clock_start)/DBLE(clock_rate)

    !create Hessenberg problem
    DO i=1,d+1
      CALL dlarnv(2, iseed, n, x)
      CALL dlagge(n, n, 1, n-1, x, p(1,n*(i-1)+1), n, iseed, work, info)
    ENDDO
    DO i=1,d+1
      ncoeff(i)=dlange('F',n,n,p(1,n*(i-1)+1),n,x)
    ENDDO
    DEALLOCATE(work, x)

    !solve problem using Hessenberg code
    CALL SYSTEM_CLOCK(count_rate=clock_rate)
    CALL SYSTEM_CLOCK(COUNT=clock_start)
    CALL dhslm(p, xr, xi, yr, yi, er, ei, berr, ncoeff, iseed, d, n)
    CALL SYSTEM_CLOCK(COUNT=clock_stop)
     
    timeStats(j,2) = DBLE(clock_stop-clock_start)/DBLE(clock_rate)
    
    DEALLOCATE(p, xr, xi, yr, yi, berr, er, ei, ncoeff, cond, ferr)

  END DO
  d = 2 * d
  !print results
  WRITE(1,'(20G15.4)', advance='no') SUM(timeStats(:,1))/m
  WRITE(1, '(A)', advance='no') ', '
  WRITE(1,'(20G15.4)', advance='no') SUM(timeStats(:,2))/m
  WRITE(1, *) 
END DO
  
 CLOSE(UNIT=1)
END PROGRAM cc_test_driver
