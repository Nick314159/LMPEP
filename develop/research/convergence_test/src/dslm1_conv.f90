SUBROUTINE dslm1_conv(p, deg, er, ei, berr, exacteigs, error)
IMPLICIT NONE
!scalar arguments
INTEGER, INTENT(IN)                 :: deg
!array arguments
DOUBLE PRECISION, INTENT(IN)        :: p(*)
DOUBLE PRECISION, INTENT(OUT)       :: berr(*), er(*), ei(*), error(*)
DOUBLE COMPLEX, DIMENSION(deg), INTENT(IN)          :: exacteigs(*)
!local scalars
LOGICAL                             :: check
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
EXTERNAL                            :: dsstart, dslcorr1, dzslcorr1
!external functions
DOUBLE PRECISION                    :: dzmod, dznrm2
EXTERNAL                            :: dzmod, dznrm2

!alpha
DO i=1,deg+1
    alpha(i)=dabs(p(i))
ENDDO
!initial estimates
CALL dsstart(alpha, deg, er, ei)
!Laguerre's Method
DO i=1,deg
    check=.False.
    DO it=1,itmax
        tol=eps*dzmod(er(i),ei(i))
        IF(dabs(ei(i))<tol) THEN
            CALL dslcorr1(p, alpha, tol, deg, i, check, er, ei, berr(i))
        ELSE
            CALL dzslcorr1(p, alpha, tol, deg, i, check, er, ei, berr(i))
        ENDIF
        IF(check) THEN
            EXIT
        ENDIF
    ENDDO
    error(it) = dznrm2(deg, dcmplx(er(1:deg), ei(1:deg))-exacteigs(1:deg), 1)/dznrm2(deg, exacteigs(1:deg), 1)
ENDDO
RETURN
END SUBROUTINE dslm1_conv
