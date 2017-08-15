PROGRAM test
  USE poly_zeroes
  IMPLICIT NONE
  INTEGER           :: clock,clock_rate,clock_start,clock_stop
  INTEGER           :: i, j, nitmax, iter
  REAL (KIND=dp), DIMENSION(:), ALLOCATABLE     ::  radius
  REAL (KIND=dp)                                ::  new_root, eps, big, small, aux, ru, ri
  COMPLEX (KIND=dp), DIMENSION(:), ALLOCATABLE  ::  root, poly
  LOGICAL, DIMENSION(:), ALLOCATABLE            ::  err

INTEGER(KIND=4)                             ::  deg, startDegree, maxDegree
CHARACTER(LEN=100)                          :: arg

  eps    = EPSILON(1.0D0)
  small  = TINY(1.0D0)
  big    = HUGE(1.0D0)
  nitmax = 60

CALL getarg(1,arg)
READ(arg, *) startDegree
CALL getarg(2,arg)
READ(arg, *) maxDegree

deg=startDegree
DO WHILE(deg<maxDegree)
  CALL SYSTEM_CLOCK(COUNT=clock)
  ALLOCATE(radius(1:deg),root(1:deg),poly(0:deg),err(deg+1))
  do i=0,deg
    call random_number(ru)
    call random_number(ri)
    poly(i)=dcmplx(0.5-ru,0.5-ri)
  enddo
  call system_clock(count_rate=clock_rate)
  call system_clock(count=clock_start)
  CALL polzeros (deg, poly, eps, big, small, nitmax, root, radius, err, iter)
  call system_clock(count=clock_stop)
  write(*,'(20G13.4)') 'PZEROS TIME =', dble(clock_stop-clock_start)/dble(clock_rate)
  deg=2*deg
  DEALLOCATE(radius,root,poly,err)
ENDDO
END PROGRAM test

