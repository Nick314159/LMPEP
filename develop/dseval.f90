!************************************************************************
!			     SUBROUTINE DSEVAL				*
!           Authors: Thomas R. Cameron, Nikolas I. Steckley             *
!************************************************************************
! Evaluates scalar polynomial with real coefficients at real number.    *
!************************************************************************
! Input Variables:                                                      *
!   p: REAL(8) array of dimension (d+1),                                *
!       contains polynomial coefficients from constant to leading.      *
!   t: REAL(8), number that polynomial is evaluated at.                 *
!   deg: INTEGER(4) degree of the polynomial                            *
!   der: INTEGER(1), derivative to be taken (0,1,2)                     *
!                                                                       *
! Output Variables:                                                     *
!   a: REAL(dp), return value.                                          *
!                                                                       *
! MEMORY: O(deg), FLOPS: O(deg)                                         *
!************************************************************************
SUBROUTINE dseval(p, t, deg, der, a)
IMPLICIT NONE
!scalar arguments
INTEGER(1), INTENT(IN)      :: der
INTEGER(4), INTENT(IN)      :: deg
REAL(8),    INTENT(IN)      :: t
REAL(8),    INTENT(INOUT)   :: a
!array arguments
REAL(8),    INTENT(IN)      :: p(*)
!local scalars
INTEGER :: k

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
