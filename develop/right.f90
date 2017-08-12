!>\brief <b>  Given as input the integer I and the vector H of logical, compute the the minimum integer IR such that IR>I and H(IL) is TRUE. </b>
!>\note Borrowed from \ref bini
!>\param[in] n
!>\verbatim Length of vector h \endverbatim
!>\param[in] h
!>\verbatim Vector of Logical \endverbatim
!>\param[in] i
!>\verbatim Integer \endverbatim
!>\param[out] ir
!>\verbatim Minimum integer such that ir>i, h(ir)=.TRUE. \endverbatim
SUBROUTINE right(n, h, i, ir)
IMPLICIT NONE
INTEGER, INTENT(IN)   :: n, i
INTEGER, INTENT(OUT)  :: ir
LOGICAL, INTENT(IN)             :: h(*)

DO ir = i+1, n
  IF (h(ir)) RETURN
ENDDO
RETURN
END SUBROUTINE right

