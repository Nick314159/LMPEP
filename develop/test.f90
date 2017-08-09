PROGRAM test
USE util
INTEGER(KIND=in4)                               :: deg, it, itmax, startDegree, maxDegree
INTEGER(KIND=in4)                               :: clock, clock_rate, clock_start, clock_stop
REAL(KIND=re8)                                  :: a, t
REAL(KIND=re8), DIMENSION(:), ALLOCATABLE       :: er, ei, p
COMPLEX(KIND=re8)                               :: ac, tc
COMPLEX(KIND=re8), DIMENSION(:), ALLOCATABLE    :: pc
CHARACTER(LEN=10)                               :: dt
CHARACTER(LEN=100)                              :: arg

CALL GETARG(1,arg)
READ(arg, *) startDegree
CALL GETARG(2,arg)
READ(arg, *) maxDegree
CALL GETARG(3,arg)
READ(arg,*) dt

itmax=100
deg=startDegree

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
    WRITE(*,*) '        DEG   |   AVG ET'
    WRITE(*,*) '        ----------------'
    WRITE(*,*) deg,' ', (DBLE(clock_stop-clock_start)/DBLE(clock_rate))/itmax
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
    WRITE(*,*) '        DEG   |   AVG ET'
    WRITE(*,*) '        ----------------'
    WRITE(*,*) deg,' ', (DBLE(clock_stop-clock_start)/DBLE(clock_rate))/itmax
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
    WRITE(*,*) '        DEG   |   AVG ET'
    WRITE(*,*) '        ----------------'
    WRITE(*,*) deg,' ', (DBLE(clock_stop-clock_start)/DBLE(clock_rate))/itmax
    DEALLOCATE(er,ei,p)
    deg=2*deg
ENDDO
ENDIF

END PROGRAM test
