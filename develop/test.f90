PROGRAM test
USE util
INTEGER(KIND=in4)                           :: deg, it, itmax, startDegree, maxDegree
INTEGER(KIND=in4)                           :: clock, clock_rate, clock_start, clock_stop
REAL(KIND=re8)                              :: a, t
REAL(KIND=re8), DIMENSION(:), ALLOCATABLE   :: p
CHARACTER(LEN=100)                          :: arg

CALL GETARG(1,arg)
READ(arg, *) startDegree
CALL GETARG(2,arg)
READ(arg, *) maxDegree

itmax=100
deg=startDegree
DO WHILE(deg<maxDegree)
    ALLOCATE(p(deg+1))
    CALL SYSTEM_CLOCK(count_rate=clock_rate)
    CALL SYSTEM_CLOCK(count=clock_start)
    DO it=1,itmax
        CALL RANDOM_NUMBER(p)
        CALL RANDOM_NUMBER(t)
        CALL dseval(p, t, deg, der, a)
    ENDDO
    CALL SYSTEM_CLOCK(count=clock_stop)
    WRITE(*,*) '  DEG   |   AVG ET'
    WRITE(*,*) '  ----------------'
    WRITE(*,*) deg, ' ', (DBLE(clock_stop-clock_start)/DBLE(clock_rate))/itmax
    DEALLOCATE(p)
    deg=2*deg
ENDDO

END PROGRAM test
