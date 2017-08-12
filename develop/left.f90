!>\brief <b>  Given as input the integer I and the vector H of logical, compute the the maximum integer IL such that IL<I and H(IL) is TRUE. </b>
!>\note Borrowed from \ref bini
!>\param[in] h
!>\verbatim Vector of Logical \endverbatim
!>\param[in] i
!>\verbatim Integer \endverbatim
!>\param[out] il
!>\verbatim Maximum integer such that il<i, h(il)=.TRUE. \endverbatim
SUBROUTINE left(h, i, il)
IMPLICIT NONE
INTEGER, INTENT(IN)   :: i
INTEGER, INTENT(OUT)  :: il
LOGICAL, INTENT(IN)             :: h(*)

DO il = i-1, 0, -1
  IF (h(il)) RETURN
ENDDO
RETURN
END SUBROUTINE left

