!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!   Thomas R. Cameron* and Nikolas I. Steckley#
!                                                               
!   *Dept. Mathematics and Computer Science, Davidson College
!   #Dept. DiscoverOrg LLC.
!       
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!   Last modified 15 October 2017
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!   Test Program for random polynomials of degree start to end
!   
!   The roots are computed via Aberth's method and two different
!   modified Laguerre's method. 
!   
!   Aberth's method is a modified Newton's method that results
!   in the simultaneous convergence of all roots of a polynomial,
!   iterates are implemented in a Gauss-Seidel fashion.
!   
!   The first modified Laguerre's method computes the roots one
!   at a time and deflates on previously computed roots.
!   
!   The second modified Laguerre's method computes all roots
!   simultaneously, iterates are implemented in a Gauss-Seidel
!   fashion.
!   
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!   Inpute Parameter
!   
!   1) Start Degree, default 100
!   
!   2) End Degree, default 6400
!   
!   3) Iterations, default 10
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
PROGRAM testMethods
IMPLICIT NONE
!   Input Variables
INTEGER :: startDegree, endDegree, iter, flag
CHARACTER(len=32) :: arg

!   Compute Variables
INTEGER :: clock, clock_rate, clock_start, clock_stop, deg, it
DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE :: berr, er, ei, p
DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE :: berrMethods, timeMethods

!   External Subroutines
EXTERNAL :: daruv, dslm, dslm1, dsam

!   Set Start and End Degree, and Iterations
    CALL GET_COMMAND_ARGUMENT(1,arg,status=flag)
    IF(flag==0) THEN
        READ(arg, '(I10)') startDegree
    ELSE
        startDegree=100
    ENDIF
    CALL GET_COMMAND_ARGUMENT(2,arg,status=flag)
    IF(flag==0) THEN
        READ(arg, '(I10)') endDegree
    ELSE
        endDegree=6400
    ENDIF
    CALL GET_COMMAND_ARGUMENT(3,arg,status=flag)
    IF(flag==0) THEN
        READ(arg, '(I10)') iter
    ELSE
        iter=10
    ENDIF

!   Open Results File
OPEN(UNIT=1,FILE="results/testMethods.csv")
WRITE(1,'(A)') 'Degree, DSLM Time, DSLM Berr, DSLM1 Time, DSLM1 Berr, DSAM Time, DSAM Berr'
!   Allocate Backward Error and Elapsed Time storage arrays
ALLOCATE(berrMethods(iter,3), timeMethods(iter,3))
!   Main Loop
deg=startDegree
DO WHILE(deg<=endDegree)
    WRITE(1,'(I10)',advance='no') deg
    WRITE(1,'(A)',advance='no') ','
    
    DO it=1,iter
        !   Random polynomial of degree deg
        ALLOCATE(p(deg+1))
        CALL daruv(deg+1,p)
        !   DSLM: DSLM_Simul
        ALLOCATE(berr(deg), er(deg), ei(deg))
        CALL SYSTEM_CLOCK(count_rate=clock_rate)
        CALL SYSTEM_CLOCK(count=clock_start)
        CALL dslm(p, deg, er, ei, berr)
        CALL SYSTEM_CLOCK(count=clock_stop)
        timeMethods(it,1)=DBLE(clock_stop-clock_start)/DBLE(clock_rate)
        berrMethods(it,1)=MAXVAL(berr)
        DEALLOCATE(berr, er, ei)
        !   DSLM1: DSLM_Seque
        ALLOCATE(berr(deg), er(deg), ei(deg))
        CALL SYSTEM_CLOCK(count_rate=clock_rate)
        CALL SYSTEM_CLOCK(count=clock_start)
        CALL dslm1(p, deg, er, ei, berr)
        CALL SYSTEM_CLOCK(count=clock_stop)
        timeMethods(it,2)=DBLE(clock_stop-clock_start)/DBLE(clock_rate)
        berrMethods(it,2)=MAXVAL(berr)
        DEALLOCATE(berr, er, ei)
        !   DSAM
        ALLOCATE(berr(deg), er(deg), ei(deg))
        CALL SYSTEM_CLOCK(count_rate=clock_rate)
        CALL SYSTEM_CLOCK(count=clock_start)
        CALL dsam(p, deg, er, ei, berr)
        CALL SYSTEM_CLOCK(count=clock_stop)
        timeMethods(it,3)=DBLE(clock_stop-clock_start)/DBLE(clock_rate)
        berrMethods(it,3)=MAXVAL(berr)
        DEALLOCATE(berr, er, ei)
        DEALLOCATE(p)
    ENDDO

    WRITE(1,'(ES15.2)',advance='no') sum(timeMethods(:,1))/iter
    WRITE(1,'(A)',advance='no') ','
    WRITE(1,'(ES15.2)',advance='no') sum(berrMethods(:,1))/iter
    WRITE(1,'(A)',advance='no') ','
    WRITE(1,'(ES15.2)',advance='no') sum(timeMethods(:,2))/iter
    WRITE(1,'(A)',advance='no') ','
    WRITE(1,'(ES15.2)',advance='no') sum(berrMethods(:,2))/iter
    WRITE(1,'(A)',advance='no') ','
    WRITE(1,'(ES15.2)',advance='no') sum(timeMethods(:,3))/iter
    WRITE(1,'(A)',advance='no') ','
    WRITE(1,'(ES15.2)',advance='no') sum(berrMethods(:,3))/iter
    WRITE(1,'(A)')
    
    deg=2*deg
ENDDO
!   Deallocate Backward Error and Elapsed Time storage arrays
DEALLOCATE(berrMethods, timeMethods)
!   Close Results File
CLOSE(1)
!   Call testMethods.py to display results
CALL EXECUTE_COMMAND_LINE('python testMethods.py')

END PROGRAM testMethods
