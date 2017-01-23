MODULE qepinit

	IMPLICIT NONE

!	INTEGER, PARAMETER				:: dp = SELECTED_REAL_KIND(15, 60)
!	REAL(dp), PARAMETER				:: eps = 10*EPSILON(1.d0)
	REAL(dp), PARAMETER				:: epsloose = 100*eps
	REAL(dp), PARAMETER				:: epsrelaxed = 1.d-06
!	REAL(dp), PARAMETER				:: zero = 0.0D0, one = 1.0D0

END MODULE qepinit
