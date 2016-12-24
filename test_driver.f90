PROGRAM driver
    USE poly_zeros
    USE dslmpep_subroutines
    IMPLICIT NONE

    !=======VARIABLES=======
    !Common
    !CHARACTER(*), PARAMETER :: resultsDir="/home/thomas/Documents/FORTRAN/Nick/LMPEPtests/results/"
    CHARACTER(*), PARAMETER :: resultsDir="/home/nsteckley/Documents/Personal/Cameron/LMPEP/tests/results/"
    INTEGER, DIMENSION(4) :: iseed
    REAL (KIND=dp), DIMENSION(:), ALLOCATABLE  ::  testProblem
    INTEGER :: clock, clock_rate, clock_start, clock_stop
    REAL (KIND=dp), DIMENSION(10,2) :: timingStatistics, radStats
    INTEGER :: degree
    INTEGER :: startDegree, maxDegree, jumpFactor
    CHARACTER *100 BUFFER

    !DSLMPEP ----------
    INTEGER :: i
    REAL(dp), DIMENSION(:), ALLOCATABLE :: backwardError, realRoots, imaginaryRoots
    !------------------
	
    !PZEROS -----------
    INTEGER :: iterations_p
    REAL (KIND=dp), DIMENSION(:), ALLOCATABLE :: radius
    COMPLEX (KIND=dp), DIMENSION(:), ALLOCATABLE :: root
    LOGICAL, DIMENSION(:), ALLOCATABLE :: err
    !------------------

    CALL GETARG(1,BUFFER)
    READ(BUFFER, *) startDegree
    CALL GETARG(2,BUFFER)
    READ(BUFFER, *) maxDegree
    CALL GETARG(3,BUFFER)
    READ(BUFFER, *) jumpFactor
	
    !Create iseed
    CALL SYSTEM_CLOCK(COUNT=clock)
    CALL srand(clock)
    DO i=1,4
        iseed(i)=MOD(irand(),4095)
    ENDDO
    IF(MOD(iseed(4),2)==0) THEN
        iseed(4)=iseed(4)+1
    ENDIF

    OPEN(UNIT=1,FILE=resultsDir//"output.csv")
    WRITE(1, '(A)',  advance='no') 'DEGREE,      '
    WRITE(1, '(A)',  advance='no') 'DSLMPEP TIME,    '
    WRITE(1, '(A)',  advance='no') 'DSLMPEP RAD,     '
    WRITE(1, '(A)',  advance='no') 'PZEROS TIME,     '
    WRITE(1, '(A)',  advance='no') 'PZEROS RAD' 
    WRITE(1, *)
    degree = startDegree
    DO WHILE (degree < maxDegree)
        WRITE(1, '(i6)', advance='no') degree
        WRITE(1, '(A)', advance='no') ', '
        CALL SYSTEM_CLOCK(count_rate=clock_rate)
        DO i=1,10
            !Generate random problem of degree i
            ALLOCATE(testProblem(degree+1))
            CALL dlarnv(2,iseed,degree+1,testProblem)
		
            !=======SOLVE=======
            !DSLMPEP ----------
            ALLOCATE(radius(degree),realRoots(degree),imaginaryRoots(degree))
            radius=zero
            CALL SYSTEM_CLOCK(count=clock_start) 
            CALL dslm(testProblem, realRoots, imaginaryRoots, radius, degree)
            CALL SYSTEM_CLOCK(count=clock_stop)
            timingStatistics(i,1) = DBLE(clock_stop-clock_start)/DBLE(clock_rate)
            radStats(i,1)=MAXVAL(radius)
            !------------------
		
            !PZEROS -----------
            ALLOCATE(root(1:degree),err(degree+1))
            CALL system_clock(count=clock_start)
            CALL polzeros (degree, DCMPLX(testProblem), eps, big, small, itmax, root, radius, err, iterations_p)
            CALL system_clock(count=clock_stop)
            timingStatistics(i,2) = DBLE(clock_stop-clock_start)/DBLE(clock_rate)
            radStats(i,2)=MAXVAL(radius)
            !------------------
		
            DEALLOCATE(testProblem)
            DEALLOCATE(radius,realRoots,imaginaryRoots)
            DEALLOCATE(root,err)
        END DO
        !=======SAVE RESULTS=======

        !DSLMPEP ----------
        WRITE(1,'(20G15.4)', advance='no') SUM(timingStatistics(:,1))/10
        WRITE(1, '(A)', advance='no') ', '	
        WRITE(1,'(20G15.4)', advance='no') SUM(radStats(:,1))/10
        WRITE(1, '(A)', advance='no') ', '
        !------------------

        !PZEROS -----------
        WRITE(1, '(20G15.4)', advance='no') SUM(timingStatistics(:,2))/10
        WRITE(1, '(A)', advance='no') ', '
        WRITE(1, '(20G15.4)', advance='no') SUM(radStats(:,2))/10
        !------------------
        WRITE(1, *)
        degree = jumpFactor * degree
    END DO
    CLOSE(UNIT=1)
END PROGRAM driver
