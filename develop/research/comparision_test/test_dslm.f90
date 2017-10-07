PROGRAM test_dslm
IMPLICIT NONE
LOGICAL                                     :: conv
INTEGER                                     :: it, itmax
INTEGER                                     :: deg, k, startDegree, maxDegree
INTEGER                                     :: clock_rate, clock_start, clock_stop
DOUBLE PRECISION                            :: a, t
DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE :: alpha, berr, berr_max, er, ei, p
DOUBLE COMPLEX                              :: ac, tc
CHARACTER(LEN=10)                           :: dt
CHARACTER(LEN=100)                          :: arg
!intrinsic subroutines
INTRINSIC                                   :: dabs, dble, dcmplx, getarg, maxval, random_number, system_clock, sum
!external subroutines
EXTERNAL                                    :: dseval, drevseval, dzseval, dzrevseval, dsstart, dslcorr, dzslcorr, dslm

CALL init_random_seed()

CALL getarg(1,arg)
READ(arg, *) startDegree
CALL getarg(2,arg)
READ(arg, *) maxDegree
CALL getarg(3,arg)
READ(arg,*) dt

itmax=10
deg=startDegree

OPEN(UNIT=1,FILE="results/results.csv")
WRITE(1,'(A)') dt

IF (dt=='DSEVAL') THEN
DO WHILE(deg<=maxDegree)
    ALLOCATE(p(deg+1))
    CALL system_clock(count_rate=clock_rate)
    CALL system_clock(count=clock_start)
    DO it=1,itmax
        CALL daruv(deg+1,p)
        CALL random_number(t)
        CALL dseval(p, t, deg, 0, a)
    ENDDO
    CALL system_clock(count=clock_stop)
    WRITE(1,'(I10)', advance='no') deg
    WRITE(1,'(A)', advance='no') ','
    WRITE(1,'(ES15.2)') (dble(clock_stop-clock_start)/dble(clock_rate))/itmax
    DEALLOCATE(p)
    deg=2*deg
ENDDO
ELSEIF(dt=='DREVSEVAL') THEN
DO WHILE(deg<=maxDegree)
    ALLOCATE(p(deg+1))
    CALL system_clock(count_rate=clock_rate)
    CALL system_clock(count=clock_start)
    DO it=1,itmax
        CALL daruv(deg+1,p)
        CALL random_number(t)
        CALL drevseval(p, t**(-1), deg, 0, a)
    ENDDO
    CALL system_clock(count=clock_stop)
    WRITE(1,'(I10)', advance='no') deg
    WRITE(1,'(A)', advance='no') ','
    WRITE(1,'(ES15.2)') (dble(clock_stop-clock_start)/dble(clock_rate))/itmax
    DEALLOCATE(p)
    deg=2*deg
ENDDO
ELSEIF(dt=='DZSEVAL') THEN
DO WHILE(deg<=maxDegree)
    ALLOCATE(p(deg+1))
    CALL system_clock(count_rate=clock_rate)
    CALL system_clock(count=clock_start)
    DO it=1,itmax
        CALL daruv(deg+1,p)
        CALL random_number(a)
        CALL random_number(t)
        tc=dcmplx(a,t)
        CALL dzseval(p, tc, deg, 0, ac)
    ENDDO
    CALL system_clock(count=clock_stop)
    WRITE(1,'(I10)', advance='no') deg
    WRITE(1,'(A)', advance='no') ','
    WRITE(1,'(ES15.2)') (dble(clock_stop-clock_start)/dble(clock_rate))/itmax
    DEALLOCATE(p)
    deg=2*deg
ENDDO
ELSEIF(dt=='DZREVSEVAL') THEN
DO WHILE(deg<=maxDegree)
    ALLOCATE(p(deg+1))
    CALL system_clock(count_rate=clock_rate)
    CALL system_clock(count=clock_start)
    DO it=1,itmax
        CALL daruv(deg+1,p)
        CALL random_number(a)
        CALL random_number(t)
        tc=dcmplx(a,t)
        CALL dzrevseval(p, tc**(-1), deg, 0, ac)
    ENDDO
    CALL system_clock(count=clock_stop)
    WRITE(1,'(I10)', advance='no') deg
    WRITE(1,'(A)', advance='no') ','
    WRITE(1,'(ES15.2)') (dble(clock_stop-clock_start)/dble(clock_rate))/itmax
    DEALLOCATE(p)
    deg=2*deg
ENDDO
ELSEIF(dt=='DSSTART') THEN
DO WHILE(deg<=maxDegree)
    ALLOCATE(alpha(deg+1),er(deg),ei(deg),p(deg+1))
    CALL system_clock(count_rate=clock_rate)
    CALL system_clock(count=clock_start)
    DO it=1,itmax
        CALL daruv(deg+1,p)
        DO k=1,deg+1
            alpha(k)=dabs(p(k))
        ENDDO
        CALL dsstart(alpha, deg, er, ei)
    ENDDO
    CALL system_clock(count=clock_stop)
    WRITE(1,'(I10)', advance='no') deg
    WRITE(1,'(A)', advance='no') ','
    WRITE(1,'(ES15.2)') (dble(clock_stop-clock_start)/dble(clock_rate))/itmax
    DEALLOCATE(alpha,er,ei,p)
    deg=2*deg
ENDDO
ELSEIF(dt=='DSLCORR') THEN
DO WHILE(deg<=maxDegree)
    ALLOCATE(alpha(deg+1),er(deg),ei(deg),p(deg+1))
    CALL system_clock(count_rate=clock_rate)
    CALL system_clock(count=clock_start)
    DO it=1,itmax
        CALL daruv(deg+1,p)
        DO k=1,deg+1
            alpha(k)=dabs(p(k))
        ENDDO
        CALL daruv(deg,er)
        CALL daruv(deg,ei)
        CALL dslcorr(p, alpha, 2*10**(-15), deg, deg, conv, er, ei, t)
    ENDDO
    CALL system_clock(count=clock_stop)
    WRITE(1,'(I10)', advance='no') deg
    WRITE(1,'(A)', advance='no') ','
    WRITE(1,'(ES15.2)') (dble(clock_stop-clock_start)/dble(clock_rate))/itmax
    DEALLOCATE(alpha,er,ei,p)
    deg=2*deg
ENDDO
ELSEIF(dt=='DZSLCORR') THEN
DO WHILE(deg<=maxDegree)
    ALLOCATE(alpha(deg+1),er(deg),ei(deg),p(deg+1))
    CALL system_clock(count_rate=clock_rate)
    CALL system_clock(count=clock_start)
    DO it=1,itmax
        CALL daruv(deg+1,p)
        DO k=1,deg+1
            alpha(k)=dabs(p(k))
        ENDDO
        CALL daruv(deg,er)
        CALL daruv(deg,ei)
        CALL dzslcorr(p, alpha, 2*10**(-15), deg, deg, conv, er, ei, t)
    ENDDO
    CALL system_clock(count=clock_stop)
    WRITE(1,'(I10)', advance='no') deg
    WRITE(1,'(A)', advance='no') ','
    WRITE(1,'(ES15.2)') (dble(clock_stop-clock_start)/dble(clock_rate))/itmax
    DEALLOCATE(alpha,er,ei,p)
    deg=2*deg
ENDDO
ELSEIF(dt=='DSLM') THEN
DO WHILE(deg<=maxDegree)
    ALLOCATE(berr(deg),berr_max(itmax),er(deg),ei(deg),p(deg+1))
    CALL system_clock(count_rate=clock_rate)
    CALL system_clock(count=clock_start)
    DO it=1,itmax
        CALL daruv(deg+1,p)
        CALL dslm(p, deg, er, ei, berr)
        berr_max(it)=maxval(berr)
    ENDDO
    CALL system_clock(count=clock_stop)
    WRITE(1,'(I10)', advance='no') deg
    WRITE(1,'(A)', advance='no') ','
    WRITE(1,'(ES15.2)', advance='no') (dble(clock_stop-clock_start)/dble(clock_rate))/itmax
    WRITE(1,'(A)', advance='no') ','
    WRITE(1,'(ES15.2)') sum(berr_max)/itmax
    DEALLOCATE(berr,berr_max,er,ei,p)
    deg=2*deg
ENDDO
ENDIF

 CLOSE(UNIT=1)

CALL EXECUTE_COMMAND_LINE('python test_dslm.py')

END PROGRAM test_dslm
