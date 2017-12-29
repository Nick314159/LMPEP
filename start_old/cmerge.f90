!>\brief <b>Given the upper convex hulls of two consecutive sets of pairs (j,A(j)), compute the upper convex hull of their union</b>
!>\note Borrowed from \ref bini
!>\param[in] n
!>\verbatim Size of the vector a \endverbatim
!>\param[in] a
!>\verbatim vector defining the points (j,A(j))  \endverbatim
!>\param[in] i
!>\verbatim abscissa of the common vertex of the two sets  \endverbatim
!>\param[in] m
!>\verbatim  the number of elements of each set is M+1 \endverbatim
!>\param[in,out] h
!>\verbatim vector defining the vertices of the convex hull, i.e.,H(j) is .TRUE. if (j,A(j)) is a vertex of the convex hull. This vector is used also as output.\endverbatim
!************************************************************************
SUBROUTINE cmerge(n, a, i, m, h)
IMPLICIT NONE
INTEGER, INTENT(IN)             :: n, m, i
LOGICAL, INTENT(INOUT)          :: h(*)
DOUBLE PRECISION, INTENT(IN)    :: a(*)

! Local variables
INTEGER                         :: ir, il, irr, ill
LOGICAL                         :: tstl, tstr
!intrinsic functions            
INTRINSIC                       :: min
!external subroutines
EXTERNAL                        :: left, right
!external functions
LOGICAL                         :: ctest
EXTERNAL                        :: ctest

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
  h(i) = .FALSE.
  DO
    IF (il == (i-m)) THEN
      tstl = .TRUE.
    ELSE
      CALL left(h, il, ill)
      tstl = ctest(a, ill, il, ir)
    ENDIF
    IF (ir == min(n, i+m)) THEN
      tstr = .TRUE.
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
