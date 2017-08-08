MODULE util

IMPLICIT NONE
!parameters
INTEGER, PARAMETER      :: in1 = SELECTED_INT_KIND(1)
INTEGER, PARAMETER      :: in4 = SELECTED_INT_KIND(6)     
INTEGER, PARAMETER      :: re8 = SELECTED_REAL_KIND(15, 307)

CONTAINS

SUBROUTINE drnum(a)
!scalar inputs
REAL(KIND=re8)  :: a
!local scalars
REAL(KIND=re8)  :: a1, a2

CALL RANDOM_NUMBER(a1)
CALL RANDOM_NUMBER(a2)
a=-2*a1+a2
END SUBROUTINE drnum

SUBROUTINE zrnum(a)
!scalar inputs
COMPLEX(KIND=re8)   :: a
!local scalars
REAL(KIND=re8)      :: a1, a2

CALL drnum(a1)
CALL drnum(a2)
a=DCMPLX(a1,a2)
END SUBROUTINE zrnum

SUBROUTINE drarr(p, n)
IMPLICIT NONE
!scalar inputs
INTEGER(KIND=in4), INTENT(IN)   :: n
!array inputs
REAL(KIND=re8), INTENT(INOUT)   :: p(*)
!local arrays
REAL(KIND=re8), DIMENSION(n)    :: p1, p2

CALL RANDOM_NUMBER(p1)
CALL RANDOM_NUMBER(p2)
p(1:n)=-2*p1+p2
RETURN
END SUBROUTINE drarr

SUBROUTINE zrarr(p, n)
IMPLICIT NONE
!scalar inputs
INTEGER(KIND=in4), INTENT(IN)       :: n
!array inputs
COMPLEX(KIND=re8), INTENT(INOUT)    :: p(*)
!local arrays
REAL(KIND=re8), DIMENSION(n)        :: p1, p2

CALL drarr(p1,n)
CALL drarr(p2,n)
p(1:n)=DCMPLX(p1,p2)
RETURN
END SUBROUTINE zrarr

END MODULE util
