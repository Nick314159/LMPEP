PROGRAM test
INTEGER :: clock, clock_rate, clock_start, clock_stop, deg
REAL(8) :: a, t
REAL(8), DIMENSION(:), ALLOCATABLE :: p

CALL SYSTEM_CLOCK(count_rate=clock_rate)
CALL SYSTEM_CLOCK(count=clock_start)
DO deg=10,100,2
    ALLOCATE(p(deg+1))
    CALL RANDOM_NUMBER(p)
    CALL RANDOM_NUMBER(t)
    CALL dseval(p, t, a, deg, 0)
    DEALLOCATE(p)
ENDDO
CALL SYSTEM_CLOCK(count=clock_stop)
WRITE(*,*) DBLE(clock_stop-clock_start)/DBLE(clock_rate)

END PROGRAM test
