MODULE util

IMPLICIT NONE
!parameters
INTEGER, PARAMETER :: dp=KIND(0.0D0), itmax=50
REAL(dp), PARAMETER :: zero=0.0_dp, one=1.0_dp
REAL(dp), PARAMETER :: eps=EPSILON(zero), big=HUGE(zero), small=TINY(zero)
COMPLEX(dp), PARAMETER :: czero=DCMPLx(zero), cone=DCMPLX(one)
REAL(dp), PARAMETER :: epsloose = 10*eps
REAL(dp), PARAMETER :: epsrelaxed = 100*eps

CONTAINS

END MODULE util
