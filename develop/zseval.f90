
!************************************************************************
!			     SUBROUTINE ZSEVAL				*
!           Authors: Thomas R. Cameron, Nikolas I. Steckley             *
!************************************************************************
! Evaluates scalar polynomial with real coefficients at complex number. *
!************************************************************************
! Input Variables:                                                      *
!   p: REAL(8) array of dimension (d+1),                                *
!       contains polynomial coefficients from constant to leading.      *
!   t: COMPLEX(8), number that polynomial is evaluated at.              *
!   deg: INTEGER(4) degree of the polynomial                            *
!   der: INTEGER(1), derivative to be taken (0,1,2)                     *
!                                                                       *
! Output Variables:                                                     *
!   a: COMPLEX(dp), return value.                                       *
!                                                                       *
! MEMORY: O(deg), FLOPS: O(deg)                                         *
!************************************************************************
SUBROUTINE zseval(p, t, a, d, der)
IMPLICIT NONE
INTEGER, PARAMETER :: dp=KIND(0.0D0), itmax=50
!scalar arguments
INTEGER, INTENT(IN) :: d, der
COMPLEX(dp), INTENT(IN) :: t
COMPLEX(dp), INTENT(INOUT) :: a
!array arguments
REAL(dp), INTENT(IN) :: p(*)
!local scalars
INTEGER :: k

IF(der==0) THEN
  a=p(d+1)
  DO k=d,1,-1
    a=t*a+p(k)
  ENDDO
ELSEIF(der==1) THEN
  a=d*p(d+1)
  DO k=d,2,-1
    a=t*a+(k-1)*p(k)
  ENDDO
ELSE
  a=d*(d-1)*p(d+1)
  DO k=d,3,-1
    a=t*a+(k-1)*(k-2)*p(k)
  ENDDO
ENDIF
RETURN
END SUBROUTINE zseval
