PROGRAM driver
  USE poly_zeroes
  IMPLICIT NONE

  INTEGER           :: clock,clock_rate,clock_start,clock_stop
  INTEGER           :: i, j, n, nitmax, iter
  INTEGER, DIMENSION(4) :: iseed
  REAL (KIND=dp), DIMENSION(:), ALLOCATABLE     ::  radius
  REAL (KIND=dp)                                ::  new_root, eps, big, small, aux, ru, ri
  COMPLEX (KIND=dp), DIMENSION(:), ALLOCATABLE  ::  root, poly
  LOGICAL, DIMENSION(:), ALLOCATABLE            ::  err

  INTERFACE
     SUBROUTINE sort(n, x, y, e)
       USE poly_zeroes
       IMPLICIT NONE
       INTEGER, INTENT(IN)               :: n
       COMPLEX (KIND=dp), INTENT(IN OUT) :: x(:)
       REAL (KIND=dp), INTENT(IN OUT)    :: y(:)
       LOGICAL, INTENT(IN OUT)           :: e(:)
     END SUBROUTINE sort
  END INTERFACE

  eps    = EPSILON(1.0D0)
  small  = TINY(1.0D0)
  big    = HUGE(1.0D0)
  nitmax = 60

  !create iseed, used in zlarnv
  CALL SYSTEM_CLOCK(COUNT=clock)
  CALL srand(clock)
  DO i=1,4
    iseed(i)=MOD(irand(),4095)
  ENDDO
  IF(MOD(iseed(4),2)==0) THEN
    iseed(4)=iseed(4)+1
  ENDIF
  ! read degree and coefficients from the file poly.pol
  !OPEN(unit=2, file='poly.pol')
  print*, 'input degree'
  read*, n
  !READ(2,*)n
  ALLOCATE(radius(1:n),root(1:n),poly(0:n),err(n+1))
  !Do i=0,n
  !   READ(2,*) aux
  !   poly(i)=aux
  !END DO
  !do i=0,n
  !  call random_number(ru)
  !  call random_number(ri)
  !  poly(i)=dcmplx(0.5-ru,0.5-ri)
  !enddo
  call zlarnv(2,iseed,n+1,poly(0:n))
  call system_clock(count_rate=clock_rate)
  call system_clock(count=clock_start)
  CALL polzeros (n, poly, eps, big, small, nitmax, root, radius, err, iter)
  call system_clock(count=clock_stop)
  write(*,'(20G12.4)') 'PZEROS TIME =', dble(clock_stop-clock_start)/dble(clock_rate)
  CALL sort(n, root, radius, err)
  !DO i=1,n
  !   WRITE(*,*) real(root(i)), aimag(root(i))
  !END DO
  WRITE(*,*)' ITERATIONS =', iter
  WRITE(*,*)' MAX RADIUS =', maxval(radius)
END PROGRAM driver

!************************************************************************
!                             SUBROUTINE SORT                           *
!************************************************************************
!   SORT  the vector X, according to nonincreasing real parts,          *
!   the same permutation is applied to vectors Y and E.                 *
!************************************************************************
SUBROUTINE sort(n, x, y, e)
  USE poly_zeroes
  IMPLICIT NONE
  INTEGER, INTENT(IN)               :: n
  COMPLEX (KIND=dp), INTENT(IN OUT) :: x(:)
  REAL (KIND=dp), INTENT(IN OUT)    :: y(:)
  LOGICAL, INTENT(IN OUT)           :: e(:)

  ! Local variables
  INTEGER           :: k, i, imax
  COMPLEX (KIND=dp) :: temp
  REAL (KIND=dp)    :: yt, amax, rxi
  LOGICAL           :: et

  DO k = 1, n-1
     amax = REAL(x(k))
     imax = k
     DO i = k+1,n
        rxi = REAL(x(i))
        IF (amax < rxi) THEN
           amax = rxi
           imax = i
        END IF
     END DO
     temp = x(k)
     x(k) = x(imax)
     x(imax) = temp
     yt = y(k)
     et = e(k)
     y(k) = y(imax)
     y(imax) = yt
     e(k) = e(imax)
     e(imax) = et
  END DO
  RETURN
END SUBROUTINE sort

