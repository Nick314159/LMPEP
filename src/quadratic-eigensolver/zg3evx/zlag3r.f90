SUBROUTINE ZLAG3R( Conj, n, m, alpha, beta, A, ldA, B, ldB, C, ldC,            &
                   nA, nB, nC, V, ldV, R, ldR, INFO )

! Purpose
! =======

! To compute the n by m matrix R given by

!    R = A*V*D1 + B*V*D2 + C*V*D3,
!
! or
!
!   R = A^H*V*D1^H + B^H*V*D2^H + C^H*V*D3^H,
!
! where A, B and C are n by n complex matrices, V is an n by m complex matrix
! and D1, D2 and D3 are the m by m diagonal matrices given by
!
!    D1 = diag(alpha(j)^2), D2 = diag(alpha(j)*beta(j)), D3 = diag(beta(j)^2),
!
! with alpha complex, but beta real.

! Arguments
! =========

! Conj  - (input) CHARACTER(1)
!         On entry: if Conj = 'N', compute
!                      R = A*V*D1 + B*V*D2 + C*V*D3,
!                   if Conj = 'C', compute
!                      R = A^H*V*D1^H + B^H*V*D2^H + C^H*V*D2^H,
!         Constraint: Conj = 'C', or 'N'.

! n     - (input) INTEGER
!         On entry: the order of the matrices A, B and C.
!         Constraint: n >= 0.

! m     - (input) INTEGER
!         On entry: the number of columns of the matrix V.
!         Constraint: m >= 0.

! alpha - (input) COMPLEX(KIND=wp) array, dimension (m)
!         On entry: alpha(j), j = 1, 2, ..., m, contains the elements of the
!                   vector alpha.

! beta  - (input) REAL(KIND=wp) array, dimension (m)
!         On entry: beta(j), j = 1, 2, ..., m, contains the elements of the
!                   vector beta.

! A     - (input) COMPLEX(KIND=wp) array, dimension (ldA,*)
!         Note: the second dimension of the array A must be at
!               least n if nA /= zero.
!         On entry: if nA /= zero, the n by n matrix A.

! ldA   - (input) INTEGER
!         On entry: the leading dimension of the array A.
!         Constraint: ldA >= n if nA /= zero.

! B     - (input) COMPLEX(KIND=wp) array, dimension (ldB,*)
!         Note: the second dimension of the array B must be at
!               least n if nB /= zero.
!         On entry: if nB /= zero, the n by n matrix B.

! ldB   - (input) INTEGER
!         On entry: the leading dimension of the array B.
!         Constraint: ldB >= n if nB /= zero.

! C     - (input) COMPLEX(KIND=wp) array, dimension (ldC,*)
!         Note: the second dimension of the array C must be at
!               least n if nC /= zero.
!         On entry: if nC /= zero, the n by n matrix C.

! ldC   - (input) INTEGER
!         On entry: the leading dimension of the array C.
!         Constraint: ldC >= n if nC /= zero.

! nA    - (input) REAL(KIND=wp)
!         On entry: if nA is zero then the matrix A is taken to be the zero
!                   matrix.

! nB    - (input) REAL(KIND=wp)
!         On entry: if nB is zero then the matrix B is taken to be the zero
!                   matrix.

! nC    - (input) REAL(KIND=wp)
!         On entry: if nC is zero then the matrix C is taken to be the zero
!                   matrix.

! V     - (input) COMPLEX(KIND=wp) array, dimension (ldV,*)
!         Note: the second dimension of the array V must be at least m.
!         On entry: the n by m matrix V.

! ldV   - (input) INTEGER
!         On entry: the leading dimension of the array V.
!         Constraint: ldV >= n.

! R     - (output) COMPLEX(KIND=wp) array, dimension (ldR,*)
!         Note: the second dimension of the array R must be at least m.
!         On exit: the n by m matrix R.

! ldR   - (input) INTEGER
!         On entry: the leading dimension of the array R.
!         Constraint: ldR >= n.

! INFO   - (output) INTEGER 
!          On exit: INFO=0 unless the routine detects an error.
!          INFO > 0
!             Allocation of memory failed.  INFO returns the value of the
!             STAT  flag from the Fortran ALLOCATE statement of the compiler
!             with which the routine was compiled.

!      .. Implicit None Statement ..
       IMPLICIT NONE
!      .. Parameters ..
       INTEGER, PARAMETER                :: wp = KIND(0.0D0)
       COMPLEX (KIND=wp), PARAMETER      :: CONE = (1.0E+0_wp,0.0E+0_wp)
       COMPLEX (KIND=wp), PARAMETER      :: CZERO = (0.0E+0_wp,0.0E+0_wp)
       REAL (KIND=wp), PARAMETER         :: ZERO = 0.0E+0_wp
!      .. Scalar Arguments ..
       INTEGER, INTENT (IN)              :: ldA, ldB, ldC, ldR, ldV, n, m
       INTEGER, INTENT (OUT)             :: INFO
       CHARACTER (1), INTENT (IN)        :: Conj
       REAL (KIND=wp), INTENT (IN)       :: nA, nB, nC
!      .. Array Arguments ..
       COMPLEX (KIND=wp), INTENT (IN)    :: alpha(m)
       COMPLEX (KIND=wp), INTENT (INOUT) :: A(ldA,*), B(ldB, *), C(ldC,*),     &
                                            V(ldV,*)
       COMPLEX (KIND=wp), INTENT (OUT)   :: R(ldR, *)
       REAL (KIND=wp), INTENT (IN)       :: beta(m)
!      .. Local Scalars ..
       LOGICAL                           :: noconj
       INTEGER                           :: iwarn_STAT, j
!      .. Local Arrays ..
       COMPLEX (KIND=wp), ALLOCATABLE    :: TEMP(:,:)
!      .. External Subroutines ..
       EXTERNAL                             ZGEMM
!      .. Intrinsic Functions ..
       INTRINSIC                            CONJG
!      .. Executable Statements ..
       CONTINUE

! Test input arguments
INFO = 0

IF( (n == 0).OR.(m == 0) )THEN
   RETURN
END IF

IF( Conj /= 'C' .AND. Conj /= 'c' )THEN
   noconj = .TRUE.
ELSE
   noconj = .FALSE.
END IF

ALLOCATE( TEMP(n,m), STAT = iwarn_STAT )
IF( iwarn_STAT /= 0 )THEN
   INFO = iwarn_STAT
   RETURN
END IF

IF( nA /= ZERO )THEN
   IF( noconj )THEN
      DO j = 1, m
         TEMP(1:n,j) = ( alpha(j)**2 )*V(1:n,j)
      END DO
   ELSE
      DO j = 1, m
         TEMP(1:n,j) = ( CONJG( alpha(j) )**2 )*V(1:n,j)
      END DO
   END IF

   CALL ZGEMM( Conj, 'NoTrans', n, m, n, CONE, A, ldA, TEMP, n, CZERO, R, n )

ELSE
   R(1:n,1:m) = CZERO
END IF

IF( nB /= ZERO )THEN
   IF( noconj )THEN
      DO j = 1, m
         TEMP(1:n,j) = ( alpha(j)*beta(j) )*V(1:n,j)
      END DO
   ELSE
      DO j = 1, m
         TEMP(1:n,j) = ( CONJG( alpha(j) )*beta(j) )*V(1:n,j)
      END DO
   END IF

   CALL ZGEMM( Conj, 'NoTrans', n, m, n, CONE, B, ldB, TEMP, n, CONE, R, n )

END IF

IF( nC /= ZERO )THEN
   DO j = 1, m
      TEMP(1:n,j) = ( beta(j)**2 )*V(1:n,j)
   END DO

   CALL ZGEMM( Conj, 'NoTrans', n, m, n, CONE, C, ldC, TEMP, n, CONE, R, n )

END IF
DEALLOCATE( TEMP, STAT = iwarn_STAT )

RETURN

END SUBROUTINE ZLAG3R
