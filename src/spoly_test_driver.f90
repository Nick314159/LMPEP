PROGRAM spoly_test_driver
    USE poly_zeros
    USE dslmpep_subroutines
    IMPLICIT NONE

    !=======VARIABLES=======
    !Common
    !CHARACTER(*), PARAMETER :: resultsDir="/home/thomas/Documents/FORTRAN/Nick/LMPEPtests/results/"
    CHARACTER(*), PARAMETER :: resultsDir="/home/nsteckley/Documents/Personal/Cameron/LMPEP/tests/results/"

    INTEGER :: clock, clock_rate, clock_start, clock_stop, i, j, m
    INTEGER :: degree, startDegree, maxDegree, jumpFactor
    INTEGER, DIMENSION(4) :: iseed
    REAL(dp), DIMENSION(:), ALLOCATABLE :: poly, radius
    REAL(dp), DIMENSION(:,:), ALLOCATABLE :: timeStats, radStats
    COMPLEX(dp) :: a, b, t
    CHARACTER(LEN=100) :: BUFFER

    !DSLMPEP ----------
    REAL(dp), DIMENSION(:), ALLOCATABLE :: backwardError, realRoots, imaginaryRoots
    !------------------
	
    !PZEROS -----------
    INTEGER :: iter
    LOGICAL, DIMENSION(:), ALLOCATABLE :: err
    COMPLEX(dp), DIMENSION(:), ALLOCATABLE :: root
    !------------------

    !stats info
    m=10
    ALLOCATE(timeStats(m,2), radStats(m,2))

    !get information
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
        DO i=1,m
            !Generate random problem of degree i
            ALLOCATE(poly(degree+1))
            CALL dlarnv(2,iseed,degree+1,poly)
		
            !=======SOLVE=======
            !DSLMPEP ----------
            ALLOCATE(backwardError(degree),realRoots(degree),imaginaryRoots(degree))
            backwardError=zero
            CALL SYSTEM_CLOCK(count=clock_start) 
            CALL dslm(poly, realRoots, imaginaryRoots, backwardError, degree)
            CALL SYSTEM_CLOCK(count=clock_stop)
            timeStats(i,1) = DBLE(clock_stop-clock_start)/DBLE(clock_rate)
            !------------------
		
            !PZEROS -----------
            ALLOCATE(root(1:degree),err(degree+1),radius(1:degree))
            CALL system_clock(count=clock_start)
            CALL polzeros(degree, DCMPLX(poly), eps, big, small, itmax, root, radius, err, iter)
            CALL system_clock(count=clock_stop)
            timeStats(i,2) = DBLE(clock_stop-clock_start)/DBLE(clock_rate)
            !------------------

            DO j=1,degree
              !DSLMPEP
              t=DCMPLX(realRoots(j),imaginaryRoots(j))
              IF(ZABS(t)>1) THEN
                CALL zrevseval(poly, 1/t, a, degree, 0)
                CALL zrevseval(poly, 1/t, b, degree, 1)
                radius(j)=ZABS(a)/ZABS(degree*a-b/t)
              ELSE
                CALL zseval(poly, t, a, degree, 0)
                CALL zseval(poly, t, b, degree, 1)
                radius(j)=ZABS(a)/(ZABS(t)*ZABS(b))
              ENDIF
            ENDDO
            radStats(i,1)=MAXVAL(radius)
            Do j=1,degree
              !PZEROS
              t=root(j)
              IF(ZABS(t)>1) THEN
                CALL zrevseval(poly, 1/t, a, degree, 0)
                CALL zrevseval(poly, 1/t, b, degree, 1)
                radius(j)=ZABS(a)/ZABS(degree*a-b/t)
              ELSE
                CALL zseval(poly, t, a, degree, 0)
                CALL zseval(poly, t, b, degree, 1)
                radius(j)=ZABS(a)/(ZABS(t)*ZABS(b))
              ENDIF
            ENDDO
            radStats(i,2)=MAXVAL(radius)
		
            DEALLOCATE(poly)
            DEALLOCATE(realRoots,imaginaryRoots,backwardError)
            DEALLOCATE(radius,root,err)
        END DO
        !=======SAVE RESULTS=======

        !DSLMPEP ----------
        WRITE(1,'(20G15.4)', advance='no') SUM(timeStats(:,1))/m
        WRITE(1, '(A)', advance='no') ', '
        WRITE(1,'(20G15.4)', advance='no') SUM(radStats(:,1))/m
        WRITE(1, '(A)', advance='no') ', '
        !------------------

        !PZEROS -----------
        WRITE(1, '(20G15.4)', advance='no') SUM(timeStats(:,2))/m
        WRITE(1, '(A)', advance='no') ', '
        WRITE(1, '(20G15.4)', advance='no') SUM(radStats(:,2))/m
        !------------------
        WRITE(1, *)
        degree = jumpFactor * degree
    END DO
    CLOSE(UNIT=1)
END PROGRAM spoly_test_driver
