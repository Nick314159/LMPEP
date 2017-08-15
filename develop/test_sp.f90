PROGRAM test
IMPLICIT NONE
LOGICAL                                     :: conv
INTEGER(KIND=1)                             :: it, itmax
INTEGER(KIND=4)                             :: deg, k, startDegree, maxDegree
INTEGER(KIND=8)                             :: clock_rate, clock_start, clock_stop
DOUBLE PRECISION                            :: a, t
DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE :: alpha, berr, berr_max, er, ei, p
DOUBLE COMPLEX                              :: ac, tc
CHARACTER(LEN=10)                           :: dt
CHARACTER(LEN=100)                          :: arg
!intrinsic subroutines
INTRINSIC                                   :: dabs, dble, getarg, maxval, random_number, system_clock, sum
!external subroutines
EXTERNAL                                    :: dseval, drevseval, dzseval, dzrevseval, dsstart, dslcorr, dzslcorr, dslm
!external functions
DOUBLE PRECISION                            :: dzmod
EXTERNAL                                    :: dzmod

CALL init_random_seed()

CALL getarg(1,arg)
READ(arg, *) startDegree
CALL getarg(2,arg)
READ(arg, *) maxDegree
CALL getarg(3,arg)
READ(arg,*) dt

itmax=127
deg=startDegree

OPEN(UNIT=1,FILE="results.csv")
WRITE(1,'(A)') dt

IF (dt=='DSEVAL') THEN
DO WHILE(deg<maxDegree)
    ALLOCATE(p(deg+1))
    CALL system_clock(count_rate=clock_rate)
    CALL system_clock(count=clock_start)
    DO it=1,itmax
        CALL drarr(p, deg+1)
        CALL drnum(t)
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
DO WHILE(deg<maxDegree)
    ALLOCATE(p(deg+1))
    CALL system_clock(count_rate=clock_rate)
    CALL system_clock(count=clock_start)
    DO it=1,itmax
        CALL drarr(p, deg+1)
        CALL drnum(t)
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
DO WHILE(deg<maxDegree)
    ALLOCATE(p(deg+1))
    CALL system_clock(count_rate=clock_rate)
    CALL system_clock(count=clock_start)
    DO it=1,itmax
        CALL drarr(p, deg+1)
        CALL zrnum(tc)
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
DO WHILE(deg<maxDegree)
    ALLOCATE(p(deg+1))
    CALL system_clock(count_rate=clock_rate)
    CALL system_clock(count=clock_start)
    DO it=1,itmax
        CALL drarr(p, deg+1)
        CALL zrnum(tc)
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
DO WHILE(deg<maxDegree)
    ALLOCATE(er(deg),ei(deg),p(deg+1))
    CALL system_clock(count_rate=clock_rate)
    CALL system_clock(count=clock_start)
    DO it=1,itmax
        CALL random_number(p)
        CALL dsstart(p, deg, er, ei)
    ENDDO
    CALL system_clock(count=clock_stop)
    WRITE(1,'(I10)', advance='no') deg
    WRITE(1,'(A)', advance='no') ','
    WRITE(1,'(ES15.2)') (dble(clock_stop-clock_start)/dble(clock_rate))/itmax
    DEALLOCATE(er,ei,p)
    deg=2*deg
ENDDO
ELSEIF(dt=='DSLCORR') THEN
DO WHILE(deg<maxDegree)
    ALLOCATE(alpha(deg+1),er(deg),ei(deg),p(deg+1))
    CALL system_clock(count_rate=clock_rate)
    CALL system_clock(count=clock_start)
    DO it=1,itmax
        CALL drarr(p,deg+1)
        DO k=1,deg+1
            alpha(k)=dabs(p(k))
        ENDDO
        CALL drarr(er,deg)
        CALL drarr(ei,deg)
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
DO WHILE(deg<maxDegree)
    ALLOCATE(alpha(deg+1),er(deg),ei(deg),p(deg+1))
    CALL system_clock(count_rate=clock_rate)
    CALL system_clock(count=clock_start)
    DO it=1,itmax
        CALL drarr(p, deg+1)
        DO k=1,deg+1
            alpha(k)=dabs(p(k))
        ENDDO
        CALL drarr(er,deg)
        CALL drarr(ei,deg)
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
DO WHILE(deg<maxDegree)
    ALLOCATE(berr(deg),berr_max(itmax),er(deg),ei(deg),p(deg+1))
    CALL system_clock(count_rate=clock_rate)
    CALL system_clock(count=clock_start)
    DO it=1,itmax
        CALL drarr(p,deg+1)
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

CALL EXECUTE_COMMAND_LINE('python test_sp.py')

CONTAINS

SUBROUTINE init_random_seed()
INTEGER                             :: i, n, clock
INTEGER, DIMENSION(:), ALLOCATABLE  :: seed

CALL RANDOM_SEED(size = n)
ALLOCATE(seed(n))

CALL SYSTEM_CLOCK(COUNT=clock)
seed = clock + 37 * (/ (i - 1, i = 1, n) /)
CALL RANDOM_SEED(PUT = seed)

DEALLOCATE(seed)
END SUBROUTINE init_random_seed

SUBROUTINE drnum(a)
!scalar inputs
DOUBLE PRECISION    :: a
!local scalars
DOUBLE PRECISION    :: a1, a2

CALL random_number(a1)
CALL random_number(a2)
a=-a1+2*a2
END SUBROUTINE drnum

SUBROUTINE zrnum(a)
!scalar inputs
DOUBLE COMPLEX      :: a
!local scalars
DOUBLE PRECISION    :: a1, a2

CALL drnum(a1)
CALL drnum(a2)
a=dcmplx(a1,a2)
END SUBROUTINE zrnum

SUBROUTINE drarr(p, n)
IMPLICIT NONE
!scalar inputs
INTEGER, INTENT(IN)             :: n
!array inputs
DOUBLE PRECISION, INTENT(INOUT) :: p(*)
!local arrays
DOUBLE PRECISION, DIMENSION(n)  :: p1, p2

CALL random_number(p1)
CALL random_number(p2)
p(1:n)=-p1+2*p2
RETURN
END SUBROUTINE drarr

END PROGRAM test
