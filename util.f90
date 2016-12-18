MODULE util

IMPLICIT NONE
!parameters
INTEGER, PARAMETER :: dp=kind(0.0D0)
REAL(dp), PARAMETER :: eps=EPSILON(0.0_dp), big=HUGE(0.0_dp), small=TINY(0.0_dp)
REAL(dp), PARAMETER :: zero=0.0_dp, one=1.0_dp
CONTAINS

SUBROUTINE cnvex(n, a, h)
IMPLICIT NONE
INTEGER, INTENT(IN)        :: n
LOGICAL, INTENT(OUT)       :: h(:)
REAL(dp), INTENT(IN) :: a(:)

! Local variables
INTEGER :: i, j, k, m, nj, jc

h(1:n) = .true.

! compute K such that N-2 <= 2**K < N-1
k = INT(LOG(n-2.0_dp)/LOG(2.0_dp))
IF(2**(k+1) <= (n-2)) k = k+1

! For each M=1,2,4,8,...,2**K, consider the NJ pairs of consecutive
! sets made up by M+1 points having the common vertex
! (JC,A(JC)), where JC=M*(2*J+1)+1 and J=0,...,NJ,
! NJ = MAX(0, INT((N-2-M)/(M+M))).
! Compute the upper convex hull of their union by means of subroutine CMERGE
m = 1
DO i = 0, k
  nj = MAX(0, INT((n-2-m)/(m+m)))
  DO j = 0, nj
    jc = (j+j+1)*m+1
    CALL cmerge(n, a, jc, m, h)
  ENDDO
  m = m+m
ENDDO
RETURN
END SUBROUTINE cnvex


SUBROUTINE left(h, i, il)
IMPLICIT NONE
INTEGER, INTENT(IN)  :: i
INTEGER, INTENT(OUT) :: il
LOGICAL, INTENT(IN)  :: h(:)

DO il = i-1, 0, -1
  IF (h(il)) RETURN
END DO
RETURN
END SUBROUTINE left

SUBROUTINE right(n, h, i, ir)
IMPLICIT NONE
INTEGER, INTENT(IN)  :: n, i
INTEGER, INTENT(OUT) :: ir
LOGICAL, INTENT(IN)  :: h(:)

DO ir = i+1, n
  IF (h(ir)) RETURN
ENDDO
RETURN
END SUBROUTINE right


SUBROUTINE cmerge(n, a, i, m, h)
IMPLICIT NONE
INTEGER, INTENT(IN) :: n, m, i
LOGICAL, INTENT(IN OUT) :: h(:)
REAL(dp), INTENT(IN) :: a(:)

! Local variables
INTEGER :: ir, il, irr, ill
LOGICAL :: tstl, tstr

! at the left and the right of the common vertex (I,A(I)) determine
! the abscissae IL,IR, of the closest vertices of the upper convex
! hull of the left and right sets, respectively
CALL left(h, i, il)
CALL right(n, h, i, ir)

! check the convexity of the angle formed by IL,I,IR
IF (ctest(a, il, i, ir)) THEN
  RETURN
ELSE
! continue the search of a pair of vertices in the left and right
! sets which yield the upper convex hull
  h(i) = .false.
  DO
    IF (il == (i-m)) THEN
      tstl = .true.
    ELSE
      CALL left(h, il, ill)
      tstl = ctest(a, ill, il, ir)
    ENDIF
    IF (ir == MIN(n, i+m)) THEN
      tstr = .true.
    ELSE
      CALL right(n, h, ir, irr)
      tstr = ctest(a, il, ir, irr)
    ENDIF
    h(il) = tstl
    h(ir) = tstr
    IF (tstl.AND.tstr) RETURN
    IF(.NOT.tstl) il = ill
    IF(.NOT.tstr) ir = irr
  ENDDO
ENDIF
RETURN
END SUBROUTINE cmerge

FUNCTION ctest(a, il, i, ir) RESULT(OK)
IMPLICIT NONE
INTEGER, INTENT(IN) :: i, il, ir
REAL(dp), INTENT(IN) :: a(:)
LOGICAL :: OK

! Local variables
REAL(dp) :: s1, s2
REAL(dp), PARAMETER :: tol = 0.4_dp

s1 = a(i) - a(il)
s2 = a(ir) - a(i)
s1 = s1*(ir-i)
s2 = s2*(i-il)
OK = .false.
IF(s1 > (s2+tol)) OK = .true.
RETURN
END FUNCTION ctest

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
END MODULE util
