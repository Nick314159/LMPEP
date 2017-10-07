PROGRAM dslm_modifications
IMPLICIT NONE
INTEGER, PARAMETER                              :: dp = SELECTED_REAL_KIND(15, 60)
INTEGER                                         :: clock, clock_rate, clock_start, clock_stop
INTEGER                                         :: i, j, nitmax, iter
DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE   :: time, backward_error
REAL(KIND=dp)                                   :: eps, big, small, aux, ru, ri
INTEGER                                         :: it, itmax
INTEGER                                         :: deg, startDegree, maxDegree, flag
DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE     :: berr, er, ei, p, alpha
CHARACTER(LEN=100)                              :: arg
DOUBLE PRECISION, ALLOCATABLE                   :: RESIDUALS(:,:)
DOUBLE PRECISION                                :: t2
DOUBLE COMPLEX                                  :: a, t
!intrinsic subroutines
INTRINSIC                                       :: dabs, dble, dcmplx, getarg, maxval, random_number, system_clock
INTRINSIC                                       :: epsilon, tiny, huge
!external subroutines
EXTERNAL                                        :: dslm, damvw, dseval, dzseval, drevseval, dzrevseval, daruv, init_random_seed

CALL init_random_seed()

eps    = epsilon(1.0D0)
small  = tiny(1.0D0)
big    = huge(1.0D0)
nitmax = 60
FLAG = 1

CALL getarg(1,arg)
READ(arg, '(I10)') startDegree
CALL getarg(2,arg)
READ(arg, '(I10)') maxDegree

deg=startDegree
itmax = 127
OPEN(UNIT=1,FILE="results/results.csv")
WRITE(1,'(A)') 'Degree, DSLM Time, DSLM berr, DSLM1 Time, DSLM1 berr, DSAM Time, DSAM berr'
ALLOCATE(time(itmax, 4), backward_error(itmax, 4))
DO WHILE(deg<=maxDegree)
  WRITE(1,'(I10)', advance='no') deg
  WRITE(1,'(A)', advance='no') ','

  DO it = 1, itmax

    ALLOCATE(p(deg+1), alpha(deg+1))
    CALL daruv(deg+1,p)
    DO j = 1, deg
       alpha(j)=dabs(p(j))
    END DO

    !DSLM 
    ALLOCATE(er(deg), ei(deg), berr(deg))
    CALL system_clock(count_rate=clock_rate)
    CALL system_clock(count=clock_start)
    CALL dslm(p, deg, er, ei, berr)
    CALL system_clock(count=clock_stop)
    time(it, 1)=(dble(clock_stop-clock_start)/dble(clock_rate))
    backward_error(it, 1) = maxval(berr)
    DEALLOCATE(er, ei, berr)

    !DSLM1 
    ALLOCATE(er(deg), ei(deg), berr(deg))
    CALL system_clock(count_rate=clock_rate)
    CALL system_clock(count=clock_start)
    CALL dslm1(p, deg, er, ei, berr)
    CALL system_clock(count=clock_stop)
    time(it, 2)=(dble(clock_stop-clock_start)/dble(clock_rate))
    backward_error(it, 2) = maxval(berr)
    DEALLOCATE(er, ei, berr)
   
    !DSAM 
    ALLOCATE(er(deg), ei(deg), berr(deg))
    CALL system_clock(count_rate=clock_rate)
    CALL system_clock(count=clock_start)
    CALL dsam(p, deg, er, ei, berr)
    CALL system_clock(count=clock_stop)
    time(it, 3)=(dble(clock_stop-clock_start)/dble(clock_rate))
    backward_error(it, 3) = maxval(berr)
    DEALLOCATE(er, ei, berr)

    DEALLOCATE(p, alpha)

   
  ENDDO
  
  WRITE(1,'(ES15.2)', advance='no') sum(time(:,1))/itmax
  WRITE(1,'(A)', advance='no') ','
  WRITE(1,'(ES15.2)', advance='no') sum(backward_error(:,1))/itmax
  WRITE(1,'(A)', advance='no') ','
  WRITE(1,'(ES15.2)', advance='no') sum(time(:,2))/itmax
  WRITE(1,'(A)', advance='no') ','
  WRITE(1,'(ES15.2)', advance='no') sum(backward_error(:,2))/itmax
  WRITE(1,'(A)', advance='no') ','
  WRITE(1,'(ES15.2)', advance='no') sum(time(:,3))/itmax
  WRITE(1,'(A)', advance='no') ','
  WRITE(1,'(ES15.2)', advance='no') sum(backward_error(:,3))/itmax 
  WRITE(1,'(A)')
  deg=2*deg
ENDDO
DEALLOCATE(time, backward_error)

CALL EXECUTE_COMMAND_LINE('python time_and_error.py')

END PROGRAM dslm_modifications

