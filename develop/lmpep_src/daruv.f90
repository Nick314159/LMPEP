!>\author Thomas R. Cameron, Davidson College
!>\author Nikolas I. Steckley, Portland State University
!>\date 2017
!>\brief <b> DARUV creates a random array of numbers </b>
!>\par Purpose:
!>\verbatim
!> DARUV creates N random double precious numbers and writes them to the array X
!>\endverbatim
!>\param[in] n
!>\verbatim Integer number of random numbers desired \endverbatim
!>\param[out] x
!>\verbatim Double precision array to write random numbers to\endverbatim
!***********************************************************************
SUBROUTINE daruv(n,x)
! scalar arguments
INTEGER, INTENT(IN)             :: n
! array arguments
DOUBLE PRECISION, INTENT(OUT)   :: x(*)
! local scalars
INTEGER                         :: k
DOUBLE PRECISION                :: r
! intrinsic functions
INTRINSIC                       :: random_number

DO k=1,n
    CALL random_number(r)
    x(k) = -1+2*r
ENDDO
RETURN
END SUBROUTINE
