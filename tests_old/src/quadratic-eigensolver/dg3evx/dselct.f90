    FUNCTION DSELCT(ar,ai,b)

!     Dummy SELCTG for use with DGGES and DGGESX

!     .. Implicit None Statement ..
      IMPLICIT NONE
!      .. Parameters ..
       INTEGER, PARAMETER              :: wp = KIND(0.0D0)
!     .. Function Return Value ..
      LOGICAL                          :: DSELCT
!     .. Scalar Arguments ..
      REAL (KIND=wp), INTENT (IN)      :: ai, ar, b
!     .. Executable Statements ..
      CONTINUE

      DSELCT = .FALSE.

      RETURN
    END FUNCTION DSELCT
