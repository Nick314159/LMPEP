!************************************************************************
!                             SUBROUTINE LEFT                           *
!************************************************************************
! Finds the largest integer less than the given for which the vector of *
! logicals is true.                                                     *
!************************************************************************
! Input variables:                                                      *
!     h   : vector of logical                                           *
!     i   : integer                                                     *
!                                                                       *
! Output variable:                                                      *
!     il  : maximum integer such that il<i, h(il)=.TRUE.                *
!                                                                       *
! MEMORY: O(deg), FLOPS: O(deg)                                         *
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
