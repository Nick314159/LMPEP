PROGRAM dslm_comparisons
USE poly_zeroes
IMPLICIT NONE
INTEGER                                         :: clock, clock_rate, clock_start, clock_stop
INTEGER                                         :: i, j, nitmax, iter
REAL(KIND=dp), DIMENSION(:), ALLOCATABLE        :: radius
DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE   :: time, backward_error
REAL(KIND=dp)                                   :: eps, big, small, aux, ru, ri
COMPLEX(KIND=dp), DIMENSION(:), ALLOCATABLE     :: root, poly
LOGICAL, DIMENSION(:), ALLOCATABLE              :: err
INTEGER                                         :: it, itmax
INTEGER                                         :: deg, startDegree, maxDegree, flag
DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE     :: berr, er, ei, p, p2, alpha
CHARACTER(LEN=100)                              :: arg
DOUBLE PRECISION, ALLOCATABLE                   :: REIGS(:), IEIGS(:), RESIDUALS(:,:)
INTEGER, ALLOCATABLE                            :: ITS(:)
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
itmax = 10
OPEN(UNIT=1,FILE="results.csv")
WRITE(1,'(A)') 'Degree, LMPEP Time, LMPEP berr, Pzeros Time, Pzeros berr, AMVW Time, AMVW berr'
ALLOCATE(time(itmax, 3), backward_error(itmax, 3))
DO WHILE(deg<=maxDegree)
  WRITE(1,'(I10)', advance='no') deg
  WRITE(1,'(A)', advance='no') ','

  DO it = 1, itmax
    !LMPEP
    ALLOCATE(p(deg+1), er(deg), ei(deg), berr(deg), alpha(deg+1))
    CALL daruv(deg+1,p)
    DO j = 1, deg
       alpha(j)=dabs(p(j))
    END DO
    CALL system_clock(count_rate=clock_rate)
    CALL system_clock(count=clock_start)
    CALL dslm(p, deg, er, ei, berr)
    CALL system_clock(count=clock_stop)
    time(it, 1)=(dble(clock_stop-clock_start)/dble(clock_rate))
    backward_error(it, 1) = maxval(berr)

    !Pzeros
    ALLOCATE(poly(0:deg), radius(1:deg), root(1:deg), err(deg+1)) 
    DO i= 0, deg
      poly(i)=dcmplx(p(i+1), 0.0D0)
    END DO
    CALL system_clock(count_rate=clock_rate)
    CALL system_clock(count=clock_start)
    CALL polzeros(deg, poly, eps, big, small, nitmax, root, radius, err, iter)
    CALL system_clock(count=clock_stop)
    time(it, 2)=(dble(clock_stop-clock_start)/dble(clock_rate))
    DO j = 1, deg
      t =  root(j)
      t2 = zabs(t)
      IF (t2>1) THEN
        t = t**(-1)
        t2 = t2**(-1)
        CALL dzrevseval(p, t, deg, 0, a)
        CALL drevseval(alpha, t2, deg, 0, berr(j))
      ELSE
        CALL dzseval(p, t, deg, 0, a)
        CALL dseval(alpha, t2, deg, 0, berr(j)) 
      END IF
      berr(j) = zabs(a)/berr(j)
    END DO
    backward_error(it, 2) = maxval(berr)
    DEALLOCATE(er, ei)
    DEALLOCATE(poly, radius ,root, err)

    !AMVW
    ALLOCATE(REIGS(deg),IEIGS(deg),ITS(deg), p2(deg))
    DO i=1,deg
        p2(deg-i+1)=p(i)/p(deg+1)
    ENDDO
    CALL system_clock(count_rate=clock_rate)
    CALL system_clock(count=clock_start)
    CALL damvw(deg, p2, REIGS, IEIGS, ITS, FLAG)
    CALL system_clock(count=clock_stop)  
    time(it, 3) = dble(clock_stop-clock_start)/dble(clock_rate)
    DO j = 1, deg
      t =  dcmplx(REIGS(j), IEIGS(j))
      t2 = zabs(t)
      IF (t2>1) THEN
        t = t**(-1)
        t2 = t2**(-1)
        CALL dzrevseval(p, t, deg, 0, a)
        CALL drevseval(alpha, t2, deg, 0, berr(j))
      ELSE
        CALL dzseval(p, t, deg, 0, a)
        CALL dseval(alpha, t2, deg, 0, berr(j)) 
      END IF
      berr(j) = zabs(a)/berr(j)
    END DO
    backward_error(it, 3) = maxval(berr)
    DEALLOCATE(p, p2, REIGS,IEIGS,ITS, alpha, berr)
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
  WRITE(1,'(ES15.2)') sum(backward_error(:,3))/itmax
  deg=2*deg
ENDDO
DEALLOCATE(time, backward_error)

CALL EXECUTE_COMMAND_LINE('python dslm_comparisons.py')

END PROGRAM dslm_comparisons

