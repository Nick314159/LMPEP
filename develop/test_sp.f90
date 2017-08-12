PROGRAM test
USE util
IMPLICIT NONE
LOGICAL                                     :: conv
INTEGER(KIND=1)                             :: it, itmax
INTEGER(KIND=4)                             :: deg, k, startDegree, maxDegree
INTEGER(KIND=8)                             :: clock_rate, clock_start, clock_stop
DOUBLE PRECISION                            :: a, t
DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE :: alpha, berr, er, ei, p
DOUBLE COMPLEX                              :: ac, tc
DOUBLE COMPLEX, DIMENSION(:), ALLOCATABLE   :: pc
CHARACTER(LEN=10)                           :: dt
CHARACTER(LEN=100)                          :: arg
!external functions
DOUBLE PRECISION                            :: dzmod
EXTERNAL                                    :: dzmod

CALL GETARG(1,arg)
READ(arg, *) startDegree
CALL GETARG(2,arg)
READ(arg, *) maxDegree
CALL GETARG(3,arg)
READ(arg,*) dt

itmax=127
deg=startDegree

OPEN(UNIT=1,FILE="results.csv")
WRITE(1,'(A)') '        DEG   |   AVG ET'

IF (dt=='DSEVAL') THEN
DO WHILE(deg<maxDegree)
    ALLOCATE(p(deg+1))
    CALL SYSTEM_CLOCK(count_rate=clock_rate)
    CALL SYSTEM_CLOCK(count=clock_start)
    DO it=1,itmax
        CALL drarr(p, deg+1)
        CALL drnum(t)
        CALL dseval(p, t, deg, 0, a)
    ENDDO
    CALL SYSTEM_CLOCK(count=clock_stop)
    WRITE(1,'(I10)', advance='no') deg
    WRITE(1,'(A)', advance='no') ','
    WRITE(1,'(ES15.2)') (DBLE(clock_stop-clock_start)/DBLE(clock_rate))/itmax
    DEALLOCATE(p)
    deg=2*deg
ENDDO
ELSEIF(dt=='DREVSEVAL') THEN
DO WHILE(deg<maxDegree)
    ALLOCATE(p(deg+1))
    CALL SYSTEM_CLOCK(count_rate=clock_rate)
    CALL SYSTEM_CLOCK(count=clock_start)
    DO it=1,itmax
        CALL drarr(p, deg+1)
        CALL drnum(t)
        CALL drevseval(p, t**(-1), deg, 0, a)
    ENDDO
    CALL SYSTEM_CLOCK(count=clock_stop)
    WRITE(1,'(I10)', advance='no') deg
    WRITE(1,'(A)', advance='no') ','
    WRITE(1,'(ES15.2)') (DBLE(clock_stop-clock_start)/DBLE(clock_rate))/itmax
    DEALLOCATE(p)
    deg=2*deg
ENDDO
ELSEIF(dt=='ZSEVAL') THEN
DO WHILE(deg<maxDegree)
    ALLOCATE(pc(deg+1))
    CALL SYSTEM_CLOCK(count_rate=clock_rate)
    CALL SYSTEM_CLOCK(count=clock_start)
    DO it=1,itmax
        CALL zrarr(pc, deg+1)
        CALL zrnum(tc)
        CALL zseval(pc, tc, deg, 0, ac)
    ENDDO
    CALL SYSTEM_CLOCK(count=clock_stop)
    WRITE(1,'(I10)', advance='no') deg
    WRITE(1,'(A)', advance='no') ','
    WRITE(1,'(ES15.2)') (DBLE(clock_stop-clock_start)/DBLE(clock_rate))/itmax
    DEALLOCATE(pc)
    deg=2*deg
ENDDO
ELSEIF(dt=='ZREVSEVAL') THEN
DO WHILE(deg<maxDegree)
    ALLOCATE(pc(deg+1))
    CALL SYSTEM_CLOCK(count_rate=clock_rate)
    CALL SYSTEM_CLOCK(count=clock_start)
    DO it=1,itmax
        CALL zrarr(pc, deg+1)
        CALL zrnum(tc)
        CALL zrevseval(pc, tc**(-1), deg, 0, ac)
    ENDDO
    CALL SYSTEM_CLOCK(count=clock_stop)
    WRITE(1,'(I10)', advance='no') deg
    WRITE(1,'(A)', advance='no') ','
    WRITE(1,'(ES15.2)') (DBLE(clock_stop-clock_start)/DBLE(clock_rate))/itmax
    DEALLOCATE(pc)
    deg=2*deg
ENDDO
ELSEIF(dt=='DSSTART') THEN
DO WHILE(deg<maxDegree)
    ALLOCATE(er(deg),ei(deg),p(deg+1))
    CALL SYSTEM_CLOCK(count_rate=clock_rate)
    CALL SYSTEM_CLOCK(count=clock_start)
    DO it=1,itmax
        CALL RANDOM_NUMBER(p)
        CALL dsstart(p, deg, er, ei)
    ENDDO
    CALL SYSTEM_CLOCK(count=clock_stop)
    WRITE(1,'(I10)', advance='no') deg
    WRITE(1,'(A)', advance='no') ','
    WRITE(1,'(ES15.2)') (DBLE(clock_stop-clock_start)/DBLE(clock_rate))/itmax
    DEALLOCATE(er,ei,p)
    deg=2*deg
ENDDO
ELSEIF(dt=='DSLCORR') THEN
DO WHILE(deg<maxDegree)
    ALLOCATE(alpha(deg+1),er(deg),ei(deg),p(deg+1))
    CALL SYSTEM_CLOCK(count_rate=clock_rate)
    CALL SYSTEM_CLOCK(count=clock_start)
    DO it=1,itmax
        CALL drarr(p,deg+1)
        DO k=1,deg+1
            alpha(k)=DABS(p(k))
        ENDDO
        CALL RANDOM_NUMBER(er)
        CALL RANDOM_NUMBER(ei)
        CALL dslcorr(p, alpha, 2*10**(-15), conv, deg, deg, er, ei, t)
    ENDDO
    CALL SYSTEM_CLOCK(count=clock_stop)
    WRITE(1,'(I10)', advance='no') deg
    WRITE(1,'(A)', advance='no') ','
    WRITE(1,'(ES15.2)') (DBLE(clock_stop-clock_start)/DBLE(clock_rate))/itmax
    DEALLOCATE(alpha,er,ei,p)
    deg=2*deg
ENDDO
ELSEIF(dt=='ZSLCORR') THEN
DO WHILE(deg<maxDegree)
    ALLOCATE(alpha(deg+1),er(deg),ei(deg),pc(deg+1))
    CALL SYSTEM_CLOCK(count_rate=clock_rate)
    CALL SYSTEM_CLOCK(count=clock_start)
    DO it=1,itmax
        CALL zrarr(pc, deg+1)
        DO k=1,deg+1
            alpha(k)=DZMOD(DBLE(pc(k)),DIMAG(pc(k)))
        ENDDO
        CALL RANDOM_NUMBER(er)
        CALL RANDOM_NUMBER(ei)
        CALL dslcorr(pc, alpha, 2*10**(-15), conv, deg, deg, er, ei, t)
    ENDDO
    CALL SYSTEM_CLOCK(count=clock_stop)
    WRITE(1,'(I10)', advance='no') deg
    WRITE(1,'(A)', advance='no') ','
    WRITE(1,'(ES15.2)') (DBLE(clock_stop-clock_start)/DBLE(clock_rate))/itmax
    DEALLOCATE(alpha,er,ei,pc)
    deg=2*deg
ENDDO
ELSEIF(dt=='DSLM') THEN
DO WHILE(deg<maxDegree)
    ALLOCATE(berr(deg),er(deg),ei(deg),p(deg+1))
    CALL SYSTEM_CLOCK(count_rate=clock_rate)
    CALL SYSTEM_CLOCK(count=clock_start)
    DO it=1,itmax
        CALL drarr(p,deg+1)
        CALL dslm(p, deg, er, ei, berr)
        RETURN
    ENDDO
    PRINT*, MAXVAL(berr)
    CALL SYSTEM_CLOCK(count=clock_stop)
    WRITE(1,'(I10)', advance='no') deg
    WRITE(1,'(A)', advance='no') ','
    WRITE(1,'(ES15.2)') (DBLE(clock_stop-clock_start)/DBLE(clock_rate))/itmax
    DEALLOCATE(berr,er,ei,p)
    deg=2*deg
ENDDO
ENDIF

 CLOSE(UNIT=1)

CALL EXECUTE_COMMAND_LINE('python test.py')

END PROGRAM test
