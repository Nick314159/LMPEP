PROGRAM test
USE util
INTEGER(KIND=in4)                               :: deg, it, itmax, startDegree, maxDegree
INTEGER(KIND=in4)                               :: clock, clock_rate, clock_start, clock_stop
REAL(KIND=re8)                                  :: a, t
REAL(KIND=re8), DIMENSION(:), ALLOCATABLE       :: p
COMPLEX(KIND=re8)                               :: ac, tc
COMPLEX(KIND=re8), DIMENSION(:), ALLOCATABLE    :: pc
CHARACTER(LEN=1)                                :: dt
CHARACTER(LEN=100)                              :: arg

CALL GETARG(1,arg)
READ(arg, *) startDegree
CALL GETARG(2,arg)
READ(arg, *) maxDegree
CALL GETARG(3,arg)
READ(arg,*) dt

itmax=100
deg=startDegree

IF (dt=='D') THEN
DO WHILE(deg<maxDegree)
    ALLOCATE(p(deg+1))
    CALL SYSTEM_CLOCK(count_rate=clock_rate)
    CALL SYSTEM_CLOCK(count=clock_start)
    DO it=1,itmax
        CALL drarr(p, deg+1)
        CALL drnum(t)
        CALL dseval(p, t, deg, der, a)
    ENDDO
    CALL SYSTEM_CLOCK(count=clock_stop)
    WRITE(*,*) '  DEG   |   AVG ET'
    WRITE(*,*) '  ----------------'
    WRITE(*,*) deg, ' ', (DBLE(clock_stop-clock_start)/DBLE(clock_rate))/itmax
    DEALLOCATE(p)
    deg=2*deg
ENDDO
ELSEIF(dt=='Z') THEN
DO WHILE(deg<maxDegree)
    ALLOCATE(pc(deg+1))
    CALL SYSTEM_CLOCK(count_rate=clock_rate)
    CALL SYSTEM_CLOCK(count=clock_start)
    DO it=1,itmax
        CALL zrarr(pc, deg+1)
        CALL zrnum(tc)
        CALL zseval(pc, tc, deg, der, ac)
    ENDDO
    CALL SYSTEM_CLOCK(count=clock_stop)
    WRITE(*,*) '  DEG   |   AVG ET'
    WRITE(*,*) '  ----------------'
    WRITE(*,*) deg, ' ', (DBLE(clock_stop-clock_start)/DBLE(clock_rate))/itmax
    DEALLOCATE(pc)
    deg=2*deg
ENDDO
ENDIF

END PROGRAM test
