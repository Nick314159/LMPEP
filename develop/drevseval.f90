!************************************************************************
!                           SUBROUTINE DREVSEVAL				        *
!           Authors: Thomas R. Cameron, Nikolas I. Steckley             *
!                           Date: 8/9/2017                              *
!************************************************************************
! Evaluates the reversal of a real scalar polynomial at a real number   *
!************************************************************************
! Input Variables:                                                      *
!   p: REAL(re8), array of dimension (d+1),                             *
!       contains polynomial coefficients from constant to leading.      *
!   t: REAL(re8), number that polynomial is evaluated at.               *
!   deg: INTEGER(in4), degree of the polynomial                         *
!   der: INTEGER(in1), derivative to be taken (0,1,2)                   *
!                                                                       *
! Output Variables:                                                     *
!   a: REAL(re8), return value.                                         *
!                                                                       *
!   MEMORY: O(deg), FLOPS: O(deg)                                       *
!************************************************************************
SUBROUTINE drevseval(p, t, deg, der, a)
USE util
IMPLICIT NONE
!scalar arguments
INTEGER(KIND=in1), INTENT(IN)   :: der
INTEGER(KIND=in4), INTENT(IN)   :: deg
REAL(KIND=re8), INTENT(IN)      :: t
REAL(KIND=re8), INTENT(INOUT)   :: a
!array arguments
REAL(KIND=re8), INTENT(IN)      :: p(*)
!local scalars
INTEGER(KIND=in4)               :: k

IF(der==0) THEN
  a=p(1)
  DO k=2,deg+1
    a=t*a+p(k)
  ENDDO
ELSEIF(der==1) THEN
  a=deg*p(1)
  DO k=2,deg
    a=t*a+(deg-k+1)*p(k)
  ENDDO
ELSE
  a=deg*(deg-1)*p(1)
  DO k=2,deg-1
    a=t*a+(deg-k+1)*(deg-k)*p(k)
  ENDDO
ENDIF
RETURN
END SUBROUTINE drevseval
