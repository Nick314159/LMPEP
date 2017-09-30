SUBROUTINE dslm2(p, deg, er, ei, berr, exacteigs)
IMPLICIT NONE
!scalar arguments
INTEGER, INTENT(IN)                 :: deg
!array arguments
DOUBLE PRECISION, INTENT(IN)        :: p(*)
DOUBLE PRECISION, INTENT(OUT)       :: berr(*), er(*), ei(*)
DOUBLE COMPLEX, INTENT(IN)          :: exacteigs(*)
!local scalars
LOGICAL, DIMENSION(deg)             :: check
INTEGER                             :: i, it
DOUBLE PRECISION                    :: tol
!local arrays
DOUBLE PRECISION, DIMENSION(deg+1)  :: alpha
!intrinsic procedures
INTRINSIC                           :: dabs, epsilon
!parameters
INTEGER, PARAMETER                  :: itmax=60
DOUBLE PRECISION, PARAMETER         :: eps=epsilon(0.0D0)
!external subroutines
EXTERNAL                            :: dsstart, dslcorr2, dzslcorr2
!external functions
DOUBLE PRECISION                    :: dzmod
EXTERNAL                            :: dzmod

!alpha
DO i=1,deg+1
    alpha(i)=dabs(p(i))
ENDDO
!check
DO i=1,deg
    check(i)=.TRUE.
ENDDO
!initial estimates
!CALL dsstart(alpha, deg, er, ei)
DO i = 1, deg
  er(i) = real(exacteigs(i)) + rand()/10
  ei(i) = aimag(exacteigs(i)) + rand()/10
ENDDO
!Laguerre's Method
DO it=1,itmax
    DO i=1,deg
        IF(check(i)) THEN
            tol=eps*dzmod(er(i),ei(i))
            IF(dabs(ei(i))<tol) THEN
                CALL dslcorr2(p, alpha, tol, deg, i, check(i), er, ei, berr(i))
            ELSE
                CALL dzslcorr2(p, alpha, tol, deg, i, check(i), er, ei, berr(i))
            ENDIF
        ENDIF
    ENDDO
ENDDO
RETURN
END SUBROUTINE dslm2
