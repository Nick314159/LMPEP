SUBROUTINE daruv(n,x)
!scalar arguments
INTEGER, INTENT(IN)             :: n
!array arguments
DOUBLE PRECISION, INTENT(OUT)   :: x(*)
!local scalars
INTEGER                         :: k
DOUBLE PRECISION                :: r
!intrinsic functions
INTRINSIC                       :: random_number

DO k=1,n
    CALL random_number(r)
    x(k)=-1+2*r
ENDDO

END SUBROUTINE
