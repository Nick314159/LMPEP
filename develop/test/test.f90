PROGRAM test
USE poly_zeroes
IMPLICIT NONE
INTEGER(KIND=8)                             :: clock, clock_rate, clock_start, clock_stop
INTEGER                                     :: i, j, nitmax, iter
REAL(KIND=dp), DIMENSION(:), ALLOCATABLE    :: radius
REAL(KIND=dp), DIMENSION(:,:), ALLOCATABLE  :: time
REAL(KIND=dp)                               :: eps, big, small, aux, ru, ri
COMPLEX(KIND=dp), DIMENSION(:), ALLOCATABLE :: root, poly
LOGICAL, DIMENSION(:), ALLOCATABLE          :: err
INTEGER(KIND=1)                             :: it, itmax
INTEGER(KIND=4)                             :: deg, k, startDegree, maxDegree
DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE :: berr, er, ei, p
CHARACTER(LEN=100)                          :: arg
!intrinsic subroutines
INTRINSIC                                   :: dabs, dble, dcmplx, getarg, maxval, random_number, system_clock
INTRINSIC                                   :: epsilon, tiny, huge
!external subroutines
EXTERNAL                                    :: dslm
!external functions
DOUBLE PRECISION                            :: dzmod
EXTERNAL                                    :: dzmod

CALL init_random_seed()

eps    = epsilon(1.0D0)
small  = tiny(1.0D0)
big    = huge(1.0D0)
nitmax = 60

CALL getarg(1,arg)
READ(arg, *) startDegree
CALL getarg(2,arg)
READ(arg, *) maxDegree

deg=startDegree
itmax = 127
OPEN(UNIT=1,FILE="results.csv")
WRITE(1,'(A)') 'Degree, Pzeros, LMPEP'
ALLOCATE(time(itmax, 2))
DO WHILE(deg<maxDegree)
  WRITE(1,'(I10)', advance='no') deg
  WRITE(1,'(A)', advance='no') ','

  DO it = 1, itmax
    !Pzeros
    ALLOCATE(radius(1:deg),root(1:deg),poly(0:deg),err(deg+1))
    DO i=0,deg
      CALL random_number(ru)
      CALL random_number(ri)
      poly(i)=dcmplx(-1+2*ru,-1+2*ri)
    END DO
    CALL system_clock(count_rate=clock_rate)
    CALL system_clock(count=clock_start)
    CALL polzeros (deg, poly, eps, big, small, nitmax, root, radius, err, iter)
    CALL system_clock(count=clock_stop)
    time(it, 1)=(dble(clock_stop-clock_start)/dble(clock_rate))
    DEALLOCATE(poly,radius,root,err)

    !LMPEP
    ALLOCATE(p(deg+1), er(deg), ei(deg), berr(deg))
    CALL daruv(deg+1,p)
    CALL system_clock(count_rate=clock_rate)
    CALL system_clock(count=clock_start)
    CALL dslm(p, deg, er, ei, berr)
    CALL system_clock(count=clock_stop)
    time(it, 2)=(dble(clock_stop-clock_start)/dble(clock_rate))
    DEALLOCATE(p, er, ei, berr)
  ENDDO
  
  WRITE(1,'(ES15.2)', advance='no') sum(time(:,1))/dble(itmax)
  WRITE(1,'(A)', advance='no') ','
  WRITE(1,'(ES15.2)') sum(time(:,2))/dble(itmax)

  deg=2*deg
ENDDO
DEALLOCATE(time)
END PROGRAM test

