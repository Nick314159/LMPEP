!************************************************************************
!                           SUBROUTINE DSSTART      			        *
!           Authors: Thomas R. Cameron, Nikolas I. Steckley             *
!                           Date: 8/9/2017                              *
!************************************************************************
! Compute the initial esimates of the roots of a real scalar polynomial *
! using the Newton Polygon method.                                      *
!************************************************************************
! Input Variables:                                                      *
!   alpha: REAL(re8) array, absolute value of the polynomial coeffs     *
!   deg: INTEGER(in4), degree of polynomial                             *
!                                                                       *
! Output Variables:                                                     *
!   er: REAL(re8) array, real part of the initial estimates             *
!   ei: REAL(re8) array, imaginary part of the initial estiamtes        *
!                                                                       *
! Memory: O(deg), FLOPS: O(deg)                                         *
!************************************************************************
SUBROUTINE dsstart(alpha, deg, er, ei)
USE util
IMPLICIT NONE
!scalar arguments
INTEGER, INTENT(IN)             :: deg
!array arguments
REAL(KIND=re8), INTENT(IN)      :: alpha(*)
REAL(KIND=re8), INTENT(INOUT)   :: er(*), ei(*)
!parameters
REAL(KIND=re8), PARAMETER       :: pi2 = 6.2831853071795865_re8, sigma = 0.7_re8
!local scalars
INTEGER(KIND=in4)               :: c, i, iold, j, nzeros
REAL(KIND=re8)                  :: ang, r, th
!local arrays
LOGICAL, DIMENSION(deg+1)       :: h
REAL(KIND=re8), DIMENSION(deg+1):: a
!intrinsic procedures
INTRINSIC                       :: DCOS, DEXP, DLOG, DSIN

!compute log(alpha)
DO i=1,deg+1
  IF(alpha(i)>=eps) THEN
    a(i)=DLOG(alpha(i))
  ELSE
    a(i)=-one
  ENDIF
ENDDO
!compute upper convex hull
CALL cnvex(deg+1,a,h)
!compute initial estimates
iold=1; c=0; th=pi2/deg
DO i=2,deg+1
  IF(h(i)) THEN
    nzeros=i-iold
    r=DEXP((a(iold)-a(i))/nzeros)
    ang=pi2/nzeros
    DO j=1,nzeros
      er(c+j)=r*DCOS(ang*j+th*i+sigma)
      ei(c+j)=r*DSIN(ang*j+th*i+sigma)
    ENDDO
    c=c+nzeros
    iold=i
  ENDIF
ENDDO
RETURN
END SUBROUTINE dsstart
