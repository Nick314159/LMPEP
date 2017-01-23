MODULE util

IMPLICIT NONE
!parameters
INTEGER, PARAMETER :: dp=KIND(0.0D0), itmax=50
REAL(dp), PARAMETER :: zero=0.0_dp, one=1.0_dp
REAL(dp), PARAMETER :: eps=EPSILON(zero), big=HUGE(zero), small=TINY(zero)
COMPLEX(dp), PARAMETER :: czero=DCMPLx(zero), cone=DCMPLX(one)
REAL(dp), PARAMETER :: epsloose = 100*eps
REAL(dp), PARAMETER :: epsrelaxed = 1.d-06
CONTAINS

!************************************************************************
!                             SUBROUTINE CNVEX                          *
!************************************************************************
! Compute  the upper convex hull of the set (i,a(i)), i.e., the set of  *
! vertices (i_k,a(i_k)), k=1,2,...,m, such that the points (i,a(i)) lie *
! below the straight lines passing through two consecutive vertices.    *
! The abscissae of the vertices of the convex hull equal the indices of *
! the TRUE  components of the logical output vector H.                  *
! The used method requires O(nlog n) comparisons and is based on a      *
! divide-and-conquer technique. Once the upper convex hull of two       *
! contiguous sets  (say, {(1,a(1)),(2,a(2)),...,(k,a(k))} and           *
! {(k,a(k)), (k+1,a(k+1)),...,(q,a(q))}) have been computed, then       *
! the upper convex hull of their union is provided by the subroutine    *
! CMERGE. The program starts with sets made up by two consecutive       *
! points, which trivially constitute a convex hull, then obtains sets   *
! of 3,5,9... points,  up to  arrive at the entire set.                 *
! The program uses the subroutine  CMERGE; the subroutine CMERGE uses   *
! the subroutines LEFT, RIGHT and CTEST. The latter tests the convexity *
! of the angle formed by the points (i,a(i)), (j,a(j)), (k,a(k)) in the *
! vertex (j,a(j)) up to within a given tolerance TOLER, where i<j<k.    *
!************************************************************************
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

!************************************************************************
!                             SUBROUTINE LEFT                           *
!************************************************************************
! Given as input the integer I and the vector H of logical, compute the *
! the maximum integer IL such that IL<I and H(IL) is TRUE.              *
!************************************************************************
! Input variables:                                                      *
!     H   : vector of logical                                           *
!     I   : integer                                                     *
!************************************************************************
! Output variable:                                                      *
!     IL  : maximum integer such that IL<I, H(IL)=.TRUE.                *
!************************************************************************
SUBROUTINE left(h, i, il)
IMPLICIT NONE
INTEGER, INTENT(IN)  :: i
INTEGER, INTENT(OUT) :: il
LOGICAL, INTENT(IN)  :: h(:)

DO il = i-1, 0, -1
  IF (h(il)) RETURN
ENDDO
RETURN
END SUBROUTINE left

!************************************************************************
!                             SUBROUTINE RIGHT                          *
!************************************************************************
!************************************************************************
! Given as input the integer I and the vector H of logical, compute the *
! the minimum integer IR such that IR>I and H(IL) is TRUE.              *
!************************************************************************
!************************************************************************
! Input variables:                                                      *
!     N   : length of the vector H                                      *
!     H   : vector of logical                                           *
!     I   : integer                                                     *
!************************************************************************
! Output variable:                                                      *
!     IR  : minimum integer such that IR>I, H(IR)=.TRUE.                *
!************************************************************************
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

!************************************************************************
!                             SUBROUTINE CMERGE                         *
!************************************************************************
! Given the upper convex hulls of two consecutive sets of pairs         *
! (j,A(j)), compute the upper convex hull of their union                *
!************************************************************************
! Input variables:                                                      *
!     N    : length of the vector A                                     *
!     A    : vector defining the points (j,A(j))                        *
!     I    : abscissa of the common vertex of the two sets              *
!     M    : the number of elements of each set is M+1                  *
!************************************************************************
! Input/Output variable:                                                *
!     H    : vector defining the vertices of the convex hull, i.e.,     *
!            H(j) is .TRUE. if (j,A(j)) is a vertex of the convex hull  *
!            This vector is used also as output.                        *
!************************************************************************
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

!************************************************************************
!                             FUNCTION CTEST                            *
!************************************************************************
! Test the convexity of the angle formed by (IL,A(IL)), (I,A(I)),       *
! (IR,A(IR)) at the vertex (I,A(I)), up to within the tolerance         *
! TOLER. If convexity holds then the function is set to .TRUE.,         *
! otherwise CTEST=.FALSE. The parameter TOLER is set to 0.4 by default. *
!************************************************************************
! Input variables:                                                      *
!     A       : vector of double                                        *
!     IL,I,IR : integers such that IL < I < IR                          *
!************************************************************************
! Output:                                                               *
!     .TRUE. if the angle formed by (IL,A(IL)), (I,A(I)), (IR,A(IR)) at *
!            the vertex (I,A(I)), is convex up to within the tolerance  *
!            TOLER, i.e., if                                            *
!            (A(I)-A(IL))*(IR-I)-(A(IR)-A(I))*(I-IL)>TOLER.             *
!     .FALSE.,  otherwise.                                              *
!************************************************************************
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

END MODULE util
