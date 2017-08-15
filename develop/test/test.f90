PROGRAM test
  USE poly_zeroes
  IMPLICIT NONE
  INTEGER           :: clock,clock_rate,clock_start,clock_stop
  INTEGER           :: i, j, nitmax, iter
  REAL (KIND=dp), DIMENSION(:), ALLOCATABLE     ::  radius, time
  REAL (KIND=dp)                                ::  new_root, eps, big, small, aux, ru, ri
  COMPLEX (KIND=dp), DIMENSION(:), ALLOCATABLE  ::  root, poly
  LOGICAL, DIMENSION(:), ALLOCATABLE            ::  err

LOGICAL                                     :: conv
INTEGER(KIND=1)                             :: it, itmax
INTEGER(KIND=4)                             :: deg, k, startDegree, maxDegree
DOUBLE PRECISION                            :: a, t
DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE :: alpha, berr, berr_max, er, ei, p
DOUBLE COMPLEX                              :: ac, tc
CHARACTER(LEN=10)                           :: dt
CHARACTER(LEN=100)                          :: arg
!intrinsic subroutines
INTRINSIC                                   :: dabs, dble, getarg, maxval, random_number, system_clock
!external subroutines
EXTERNAL                                    :: dseval, drevseval, dzseval, dzrevseval, dsstart, dslcorr, dzslcorr, dslm
!external functions
DOUBLE PRECISION                            :: dzmod
EXTERNAL                                    :: dzmod

  eps    = EPSILON(1.0D0)
  small  = TINY(1.0D0)
  big    = HUGE(1.0D0)
  nitmax = 60

CALL getarg(1,arg)
READ(arg, *) startDegree
CALL getarg(2,arg)
READ(arg, *) maxDegree

deg=startDegree
itmax = 10
OPEN(UNIT=1,FILE="results.csv")

DO WHILE(deg<maxDegree)
  WRITE(1,'(I10)', advance='no') deg
  WRITE(1,'(A)', advance='no') ','
  !Pzeroes
  DO it = 1, itmax
  ALLOCATE(radius(1:deg),root(1:deg),poly(0:deg),err(deg+1), time(itmax))
  do i=0,deg
    call random_number(ru)
    call random_number(ri)
    poly(i)=dcmplx(0.5-ru,0.5-ri)
  enddo
  call system_clock(count_rate=clock_rate)
  call system_clock(count=clock_start)
  CALL polzeros (deg, poly, eps, big, small, nitmax, root, radius, err, iter)
  call system_clock(count=clock_stop)
  time(it)=(dble(clock_stop-clock_start)/dble(clock_rate))
  end do
  WRITE(1,'(ES15.2)', advance='no') sum(time)/dble(itmax)
  WRITE(1,'(A)', advance='no') ','
  DEALLOCATE(radius,root,err, time)
  !LMPEP
  DO it = 1, itmax
  ALLOCATE(p(deg+1), er(deg), ei(deg), berr(deg), time(itmax)
  DO i=0, deg
    p(i+1)=poly(i)
  ENDDO
  CALL system_clock(count_rate=clock_rate)
  CALL system_clock(count=clock_start)
  CALL dslm(p, deg, er, ei, berr)
  CALL system_clock(count=clock_stop)
  time(it)=(dble(clock_stop-clock_start)/dble(clock_rate))
  END DO
  WRITE(1,'(ES15.2)') sum(time)/dble(itmax)
  DEALLOCATE(poly, p, er, ei, berr, time)

  deg=2*deg
ENDDO
END PROGRAM test

