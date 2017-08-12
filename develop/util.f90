!>\author Thomas R. Cameron, Davidson College
!>\author Nikolas I. Steckley, Portland State University
!>\date 2017
!>\brief <b> Several parameters and subroutines. </b>
!>\par Purpose:
!>\verbatim
!> Defines several parameters and subroutines that are used throughout LMPEP. 
!>\endverbatim
!************************************************************************
MODULE util
IMPLICIT NONE
!parameters
DOUBLE PRECISION, PARAMETER   :: zero=0.0D0, one=1.0D0
DOUBLE PRECISION, PARAMETER   :: eps=EPSILON(zero)
DOUBLE PRECISION, PARAMETER   :: pi2 = 6.283185307179586, sigma = 0.7D0

CONTAINS

SUBROUTINE drnum(a)
!scalar inputs
DOUBLE PRECISION    :: a
!local scalars
DOUBLE PRECISION    :: a1, a2

CALL RANDOM_NUMBER(a1)
CALL RANDOM_NUMBER(a2)
a=-a1+2*a2
END SUBROUTINE drnum

SUBROUTINE zrnum(a)
!scalar inputs
DOUBLE COMPLEX      :: a
!local scalars
DOUBLE PRECISION    :: a1, a2

CALL drnum(a1)
CALL drnum(a2)
a=DCMPLX(a1,a2)
END SUBROUTINE zrnum

SUBROUTINE drarr(p, n)
IMPLICIT NONE
!scalar inputs
INTEGER, INTENT(IN)             :: n
!array inputs
DOUBLE PRECISION, INTENT(INOUT) :: p(*)
!local arrays
DOUBLE PRECISION, DIMENSION(n)  :: p1, p2

CALL RANDOM_NUMBER(p1)
CALL RANDOM_NUMBER(p2)
p(1:n)=-p1+2*p2
RETURN
END SUBROUTINE drarr

SUBROUTINE zrarr(p, n)
IMPLICIT NONE
!scalar inputs
INTEGER, INTENT(IN)       :: n
!array inputs
DOUBLE COMPLEX, INTENT(INOUT)    :: p(*)
!local arrays
DOUBLE PRECISION, DIMENSION(n)        :: p1, p2

CALL drarr(p1,n)
CALL drarr(p2,n)
p(1:n)=DCMPLX(p1,p2)
RETURN
END SUBROUTINE zrarr

END MODULE util
