
!************************************************************************
!			            SUBROUTINE ZSEVAL				                *
!           Authors: Thomas R. Cameron, Nikolas I. Steckley             *
!                       DATE: 8/7/2017                                  *
!************************************************************************
!Evaluates scalar polynomial with complex coefficients at complex number*
!************************************************************************
! Input Variables:                                                      *
!   p: COMPLEX(re8) array of dimension (d+1),                           *
!       contains polynomial coefficients from constant to leading.      *
!   t: COMPLEX(re8), number that polynomial is evaluated at.            *
!   deg: INTEGER(in4) degree of the polynomial                          *
!   der: INTEGER(in1), derivative to be taken (0,1,2)                   *
!                                                                       *
! Output Variables:                                                     *
!   a: COMPLEX(re8), return value.                                      *
!                                                                       *
! MEMORY: O(deg), FLOPS: O(deg)                                         *
!************************************************************************
SUBROUTINE zseval(p, t, deg, der, a)
USE util
IMPLICIT NONE
!scalar arguments
INTEGER(KIND=in1), INTENT(IN)       :: der
INTEGER(KIND=in4), INTENT(IN)       :: deg
COMPLEX(KIND=re8),    INTENT(IN)    :: t
COMPLEX(KIND=re8),    INTENT(INOUT) :: a
!array arguments
COMPLEX(KIND=re8),    INTENT(IN)    :: p(*)
!local scalars
INTEGER                             :: k

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
END SUBROUTINE zseval
