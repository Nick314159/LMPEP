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
INTRINSIC                                   :: dabs, dble, dcmplx, getarg, maxval, random_number, system_clock, sum
!external subroutines
EXTERNAL                                    ::  dslm

CALL init_random_seed()

CALL getarg(1,arg)
READ(arg, *) startDegree
CALL getarg(2,arg)
READ(arg, *) maxDegree

itmax=127
deg=startDegree

OPEN(UNIT=1,FILE="results.csv")

DO WHILE(deg<maxDegree)
    ALLOCATE(berr(deg),berr_max(itmax),er(deg),ei(deg),p(deg+1))
    DO it=1,itmax
        CALL daruv(deg+1,p)
        CALL dslm(p, deg, er, ei, berr)
        berr_max(it)=maxval(berr)
    ENDDO
    WRITE(1,'(I10)', advance='no') deg
    WRITE(1,'(A)', advance='no') ','
    WRITE(1,'(ES15.2)') sum(berr_max)/itmax
    DEALLOCATE(berr,berr_max,er,ei,p)
    deg=2*deg
ENDDO

 CLOSE(UNIT=1)

END PROGRAM test
