    FUNCTION ZSELCT( a, b )

!     Dummy SELCTG for use with ZGGES and ZGGESX

!     .. Implicit None Statement ..
      IMPLICIT NONE
!      .. Parameters ..
       INTEGER, PARAMETER              :: wp = KIND(0.0D0)
!     .. Function Return Value ..
      LOGICAL                          :: ZSELCT
!     .. Scalar Arguments ..
      COMPLEX (KIND=wp), INTENT (IN)   :: a, b
!     .. Executable Statements ..
      CONTINUE

      ZSELCT = .FALSE.

      RETURN
    END FUNCTION ZSELCT
