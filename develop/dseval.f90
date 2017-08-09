!************************************************************************
!                           SUBROUTINE DSEVAL				            *
!           Authors: Thomas R. Cameron, Nikolas I. Steckley             *
!                           Date: 8/7/2017                              *
!************************************************************************
! Evaluates scalar polynomial with real coefficients at real number.    *
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
SUBROUTINE dseval(p, t, deg, der, a)
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
  a=p(deg+1)
  DO k=deg,1,-1
    a=t*a+p(k)
  ENDDO
ELSEIF(der==1) THEN
  a=deg*p(deg+1)
  DO k=deg,2,-1
    a=t*a+(k-1)*p(k)
  ENDDO
ELSE
  a=deg*(deg-1)*p(deg+1)
  DO k=deg,3,-1
    a=t*a+(k-1)*(k-2)*p(k)
  ENDDO
ENDIF
RETURN
END SUBROUTINE dseval
