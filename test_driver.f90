PROGRAM driver
        !USE poly_zeros
        USE dslmpep_subroutines
        IMPLICIT NONE
        
        !=======VARIABLES=======
        !Common
        CHARACTER(*), PARAMETER :: resultsDir="/home/nsteckley/Documents/Personal/Cameron/LMPEP/DYLMPEP/tests/results"
        INTEGER, DIMENSION(4) :: iseed
       	REAL (KIND=dp), DIMENSION(:), ALLOCATABLE  ::  testProblem
       	INTEGER :: clock, clock_rate, clock_start, clock_stop
       	REAL (KIND=dp), DIMENSION(2) :: timingStatistics
        INTEGER :: maxDegree, degree
        
        !DSLMPEP ----------
        INTEGER :: i
        INTEGER, DIMENSION(:), ALLOCATABLE :: iterations_d
	REAL(dp), DIMENSION(:), ALLOCATABLE :: backwardError, realRoots, imaginaryRoots
        !------------------
	
	!PZEROS -----------
	INTEGER           :: nitmax, iterations_p
	REAL (KIND=dp), DIMENSION(:), ALLOCATABLE     ::  radius
	COMPLEX (KIND=dp), DIMENSION(:), ALLOCATABLE  ::  root
	LOGICAL, DIMENSION(:), ALLOCATABLE            ::  err	
        !------------------
        	
	!Create iseed
	CALL SYSTEM_CLOCK(COUNT=clock)
	CALL srand(clock)
	DO i=1,4
	  iseed(i)=MOD(irand(),4095)
	ENDDO
	IF(MOD(iseed(4),2)==0) THEN
	  iseed(4)=iseed(4)+1
	ENDIF
       
	!Read maxDegree from user       
        PRINT*, 'Input highest degree to which you would like to test'
        READ *, maxDegree
        
	OPEN(UNIT=1,FILE=resultsDir//"output.txt")

        DO degree = 1, maxDegree
        
		!Generate random problem of degree i
		ALLOCATE(testProblem(degree+1))
	  	CALL dlarnv(2,iseed,degree+1,testProblem)
		
		!=======SOLVE=======
		!DSLMPEP ----------
		ALLOCATE(backwardError(degree),realRoots(degree),imaginaryRoots(degree),iterations_d(degree))

		CALL SYSTEM_CLOCK(count_rate=clock_rate)
		CALL SYSTEM_CLOCK(count=clock_start)
		CALL dslm(testProblem, realRoots, imaginaryRoots, backwardError, iterations_d, degree)
		CALL SYSTEM_CLOCK(count=clock_stop)
		timingStatistics(1) = DBLE(clock_stop-clock_start)/DBLE(clock_rate)
		!------------------
		
		!PZEROS -----------
		!ALLOCATE(radius(1:degree),root(1:degree),testProblem(0:degree),err(degree+1))
		!CALL system_clock(count_rate=clock_rate)
		!CALL system_clock(count=clock_start)
		!CALL polzeros (degree, DCMPLX(testProblem), eps, big, small, nitmax, root, radius, err, iter)
		!CALL system_clock(count=clock_stop)
		!timingStatistics(2) = DBLE(clock_stop-clock_start)/DBLE(clock_rate)
		!------------------
		
		DEALLOCATE(testProblem)
		!=======SAVE RESULTS=======

		!DSLMPEP ----------
		WRITE(*,'(20G15.4)') 'DSLMPEP TIME =', timingStatistics(1)
		WRITE(*,'(20G15.4)') 'AVG. ITER    =', DBLE(SUM(iterations_d))/(degree)
		WRITE(*,'(20G15.4)') 'MAX BE       =', MAXVAL(backwardError)
		DO i=1,degree
		  WRITE(1,*)  realRoots(i), imaginaryRoots(i), backwardError(i)
		ENDDO
		 CLOSE(UNIT=1)
		!------------------
	
		!PZEROS -----------
		!WRITE(*,'(20G12.4)') 'PZEROS TIME =', timingStatistics(2)
		!CALL sort(n, root, radius, err)
		!WRITE(*,*)' ITERATIONS =', iter
		!WRITE(*,*)' MAX RADIUS =', maxval(radius)
		!------------------
		DEALLOCATE(backwardError,realRoots,imaginaryRoots,iterations_d)
		!DEALLOCATE(radius,root,testProblem,err)
		        
        END DO
        CLOSE(UNIT=1)
END PROGRAM driver
