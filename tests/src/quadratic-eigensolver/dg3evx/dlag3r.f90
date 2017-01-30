SUBROUTINE DLAG3R( Trans, n, m, k, indx, alpha, beta,                          &
                   A, ldA, B, ldB, C, ldC, nA, nB, nC, V, ldV, R, ldR, INFO )

! Purpose
! =======

! To compute the n by k matrix R given by

!    R = A*Vr*D1 + B*Vr*D2 + C*Vr*D3,
!
! or
!
!    R = A^T*Vr*D1 + B^T*Vr*D2 + C^T*Vr*D3,
!
! where A, B and C are n by n real matrices, Vr are the k columns of an
! n by m real matrix V that represent real vectors and D1, D2 and D3 are the
! k by m diagonal matrices given by
!
!    D1 = diag(alpha(j)^2), D2 = diag(alpha(j)*beta(j)), D3 = diag(beta(j)^2),
!
! with alpha and beta real.

! This routine is called by the quadratic eigenvalue problem routine dg3evx.

! Arguments
! =========

! Trans - (input) CHARACTER(1)
!         On entry: if Trans = 'N', compute
!                      R = A*V*D1 + B*V*D2 + C*V*D3,
!                   if Trans = 'T', compute
!                      R = A^T*V*D1 + B^T*V*D2 + C^T*V*D2,
!         Constraint: Trans = 'T', or 'N'.

! n     - (input) INTEGER
!         On entry: the order of the matrices A, B and C.
!         Constraint: n >= 0.

! m     - (input) INTEGER
!         On entry: the number of columns of the matrix V.
!         Constraint: m >= k.

! k     - (input) INTEGER
!         On entry: the number of columns of V that represent the real vectors,
!                   Vr.
!         Constraint: k >= 0.

! indx  - (input) INTEGER array, dimension (k)
!         On entry: the indices of the real eigenvalues.

! alpha - (input) REAL(KIND=wp) array, dimension (m)
!         On entry: alpha(j), j = indx(1), indx(2), ..., indx(k), contains the
!                   elements of the vector alpha corresponding to the columns of
!                   Vr.

! beta  - (input) REAL(KIND=wp) array, dimension (m)
!         On entry: beta(j), j = indx(1), indx(2), ..., indx(k), contains the
!                   elements of the vector beta corresponding to the columns of
!                   Vr.

! A     - (input) REAL(KIND=wp) array, dimension (ldA,*)
!         Note: the second dimension of the array A must be at
!               least n if nA /= zero.
!         On entry: if nA /= zero, the n by n matrix A.

! ldA   - (input) INTEGER
!         On entry: the leading dimension of the array A.
!         Constraint: ldA >= n if nA /= zero.

! B     - (input) REAL(KIND=wp) array, dimension (ldB,*)
!         Note: the second dimension of the array B must be at
!               least n if nB /= zero.
!         On entry: if nB /= zero, the n by n matrix B.

! ldB   - (input) INTEGER
!         On entry: the leading dimension of the array B.
!         Constraint: ldB >= n if nB /= zero.

! C     - (input) REAL(KIND=wp) array, dimension (ldC,*)
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

! V     - (input) REAL(KIND=wp) array, dimension (ldV,*)
!         Note: the second dimension of the array V must be at least m.
!         On entry: the n by m matrix V.  The columns indx(1), indx(2), ...,
!                   indx(k) of V must give the matrix Vr.

! ldV   - (input) INTEGER
!         On entry: the leading dimension of the array V.
!         Constraint: ldV >= n.

! R     - (output) REAL(KIND=wp) array, dimension (ldR,*)
!         Note: the second dimension of the array R must be at least k.
!         On exit: the n by k matrix R.

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
       INTEGER, PARAMETER             :: wp = KIND(0.0D0)
       REAL (KIND=wp), PARAMETER      :: ONE = 1.0E+0_wp
       REAL (KIND=wp), PARAMETER      :: ZERO = 0.0E+0_wp
!      .. Scalar Arguments ..
       INTEGER, INTENT (IN)           :: k, ldA, ldB, ldC, ldR, ldV, n, m
       INTEGER, INTENT (OUT)          :: INFO
       CHARACTER (1), INTENT (IN)     :: Trans
       REAL (KIND=wp), INTENT (IN)    :: nA, nB, nC
!      .. Array Arguments ..
       INTEGER, INTENT (IN)           :: indx(k)
       REAL (KIND=wp), INTENT (IN)    :: alpha(m), beta(m)
       REAL (KIND=wp), INTENT (INOUT) :: A(ldA,*), B(ldB, *), C(ldC,*),        &
                                         V(ldV,*)
       REAL (KIND=wp), INTENT (OUT)   :: R(ldR, *)
!      .. Local Scalars ..
       INTEGER                        :: i, iwarn_STAT, j
!      .. Local Arrays ..
       REAL (KIND=wp), ALLOCATABLE    :: TEMP(:,:)
!      .. External Subroutines ..
       EXTERNAL                          DGEMM
!      .. Executable Statements ..
       CONTINUE

! Test input arguments
INFO = 0

IF( (n == 0).OR.(k == 0) )THEN
   RETURN
END IF

ALLOCATE( TEMP(n,k), STAT = iwarn_STAT )
IF( iwarn_STAT /= 0 )THEN
   INFO = iwarn_STAT
   RETURN
END IF

IF( nA /= ZERO )THEN
   DO j = 1, k
      i = indx(j)
      TEMP(1:n,j) = ( alpha(i)**2 )*V(1:n,i)
   END DO

   CALL DGEMM( Trans, 'NoTrans', n, k, n, ONE, A, ldA, TEMP, n, ZERO, R, n )

ELSE
   R(1:n,1:m) = ZERO
END IF

IF( nB /= ZERO )THEN
   DO j = 1, k
      i = indx(j)
      TEMP(1:n,j) = ( alpha(i)*beta(i) )*V(1:n,i)
   END DO

   CALL DGEMM( Trans, 'NoTrans', n, k, n, ONE, B, ldB, TEMP, n, ONE, R, n )

END IF

IF( nC /= ZERO )THEN
   DO j = 1, k
      i = indx(j)
      TEMP(1:n,j) = ( beta(i)**2 )*V(1:n,i)
   END DO

   CALL DGEMM( Trans, 'NoTrans', n, k, n, ONE, C, ldC, TEMP, n, ONE, R, n )

END IF
DEALLOCATE( TEMP, STAT = iwarn_STAT )

RETURN

END SUBROUTINE DLAG3R
