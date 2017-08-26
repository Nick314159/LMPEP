SUBROUTINE DLAG3C( Trans, n, m, k, indx, alphar, alphai, beta,                  &
                   A, ldA, B, ldB, C, ldC, nA, nB, nC, V, ldV, R, ldR, INFO )

! Purpose
! =======

! To compute the n by k matrix R given by

!    R = A*Vc*D1 + B*Vc*D2 + C*Vc*D3,
!
! or
!
!    R = A^T*Vc*D1^H + B^T*Vc*D2^H + C^T*Vc*D3^H
!
! where A, B and C are n by n real matrices, Vc are the 2k columns of an
! n by m real matrix V that represent complex vectors and D1, D2 and D3 are the
! k by m diagonal matrices given by
!
!    D1 = diag(alpha(j)^2), D2 = diag(alpha(j)*beta(j)), D3 = diag(beta(j)^2),
!
! with alpha = alphar + i*alphai, where alphar, alphai and beta are real.

! This routine is called by the quadratic eigenvalue problem routine dg3evx.

! Arguments
! =========

! Trans  - (input) CHARACTER(1)
!          On entry: if Trans = 'N', compute
!                       R = A*Vc*D1 + B*Vc*D2 + C*Vc*D3,
!                    if Trans = 'T', compute
!                       R = A^T*Vc*D1^H + B^T*Vc*D2^H + C^T*Vc*D2^H,
!          Constraint: Trans = 'T', or 'N'.

! n      - (input) INTEGER
!          On entry: the order of the matrices A, B and C.
!          Constraint: n >= 0.

! m      - (input) INTEGER
!          On entry: the number of columns of the matrix V.
!          Constraint: m >= k.

! k      - (input) INTEGER
!          On entry: the number of pairs columns of V that represent the complex
!          vectors, Vc.  A column of Vc is given in the form V(j) + i*V(j+1).
!          Constraint: k >= 0.

! indx   - (input) INTEGER array, dimension (k)
!          On entry: the indices of the complex eigenvalues.

! alphar - (input) REAL(KIND=wp) array, dimension (m)
!          On entry: alphar(j), j = indx(1), indx(2), ..., indx(k), contains the
!                    real parts of the complex vector alpha corresponding to the
!                    columns of Vc.

! alphai - (input) REAL(KIND=wp) array, dimension (m)
!          On entry: alphai(j), j = indx(1), indx(2), ..., indx(k), contains the
!                    imaginary parts of the complex vector alpha corresponding
!                    to the columns of Vc.

! beta   - (input) REAL(KIND=wp) array, dimension (m)
!          On entry: beta(j), j = indx(1), indx(2), ..., indx(k), contains the
!                    elements of the vector beta corresponding to the columns of
!                    Vc.

! A      - (input) REAL(KIND=wp) array, dimension (ldA,*)
!          Note: the second dimension of the array A must be at
!                least n if nA /= zero.
!          On entry: if nA /= zero, the n by n matrix A.

! ldA    - (input) INTEGER
!          On entry: the leading dimension of the array A.
!          Constraint: ldA >= n if nA /= zero.

! B      - (input) REAL(KIND=wp) array, dimension (ldB,*)
!          Note: the second dimension of the array B must be at
!                least n if nB /= zero.
!          On entry: if nB /= zero, the n by n matrix B.

! ldB    - (input) INTEGER
!          On entry: the leading dimension of the array B.
!          Constraint: ldB >= n if nB /= zero.

! C      - (input) REAL(KIND=wp) array, dimension (ldC,*)
!          Note: the second dimension of the array C must be at
!                least n if nC /= zero.
!          On entry: if nC /= zero, the n by n matrix C.

! ldC    - (input) INTEGER
!          On entry: the leading dimension of the array C.
!          Constraint: ldC >= n if nC /= zero.

! nA     - (input) REAL(KIND=wp)
!          On entry: if nA is zero then the matrix A is taken to be the zero
!                    matrix.

! nB     - (input) REAL(KIND=wp)
!          On entry: if nB is zero then the matrix B is taken to be the zero
!                    matrix.

! nC     - (input) REAL(KIND=wp)
!          On entry: if nC is zero then the matrix C is taken to be the zero
!                    matrix.

! V      - (input) REAL(KIND=wp) array, dimension (ldV,*)
!          Note: the second dimension of the array V must be at least m.
!          On entry: the n by m matrix V.  The columns indx(1), indx(1)+1,
!          indx(2), indx(2)+1, ..., indx(k), indx(k)+1 of V must give the
!          matrix Vc.

! ldV    - (input) INTEGER
!          On entry: the leading dimension of the array V.
!          Constraint: ldV >= n.

! R      - (output) COMPLEX(KIND=wp) array, dimension (ldR,*)
!          Note: the second dimension of the array R must be at least k.
!          On exit: the n by k matrix R.
!          Note: the ith column of R corresponds to the indx(i)th complex
!                eigenvalue alpha(indx(i)).  The conjugate of this column is
!                not computed, so R contains k columns rather than 2k columns.

! ldR    - (input) INTEGER
!          On entry: the leading dimension of the array R.
!          Constraint: ldR >= n.

! INFO    - (output) INTEGER 
!           On exit: INFO=0 unless the routine detects an error.
!           INFO > 0
!              Allocation of memory failed.  INFO returns the value of the
!              STAT  flag from the Fortran ALLOCATE statement of the compiler
!              with which the routine was compiled.

!      .. Implicit None Statement ..
       IMPLICIT NONE
!      .. Parameters ..
       INTEGER, PARAMETER              :: wp = KIND(0.0D0)
       REAL (KIND=wp), PARAMETER       :: ZERO = 0.0E+0_wp
       COMPLEX (KIND=wp), PARAMETER    :: CONE = (1.0E+0_wp,0.0E+0_wp)
       COMPLEX (KIND=wp), PARAMETER    :: CZERO = (0.0E+0_wp,0.0E+0_wp)
!      .. Scalar Arguments ..
       INTEGER, INTENT (IN)            :: k, ldA, ldB, ldC, ldR, ldV, n, m
       INTEGER, INTENT (OUT)           :: INFO
       CHARACTER (1), INTENT (IN)      :: Trans
       REAL (KIND=wp), INTENT (IN)     :: nA, nB, nC
!      .. Array Arguments ..
       INTEGER, INTENT (IN)            :: indx(k)
       REAL (KIND=wp), INTENT (IN)     :: alphai(m), alphar(m), beta(m)
       REAL (KIND=wp), INTENT (INOUT)  :: A(ldA,*), B(ldB, *), C(ldC,*),        &
                                          V(ldV,*)
       COMPLEX (KIND=wp), INTENT (OUT) :: R(ldR, *)
!      .. Local Scalars ..
       INTEGER                         :: i, iwarn_STAT, j
       CHARACTER (1)                   :: Conj
!      .. Local Arrays ..
       COMPLEX (KIND=wp), ALLOCATABLE  :: alpha(:), CTEMP(:,:), TEMP(:,:)
!      .. External Subroutines ..
       EXTERNAL                           ZGEMM
!      .. Intrinsic Functions ..
       INTRINSIC                          CMPLX, CONJG
!      .. Executable Statements ..
       CONTINUE

! Test input arguments
INFO = 0

IF( (n == 0).OR.(k == 0) )THEN
   RETURN
END IF

ALLOCATE( alpha(k), CTEMP(n,k), TEMP(n,n), STAT = iwarn_STAT )
IF( iwarn_STAT /= 0 )THEN
   INFO = iwarn_STAT
   RETURN
END IF

IF( Trans == 'N' .OR. TRANS == 'n' )THEN
   Conj = 'N'
ELSE
   Conj = 'C'
END IF

DO j = 1, k
   i = indx(j)
   IF( Conj == 'N' )THEN
      alpha(j) = CMPLX( alphar(i), alphai(i), KIND = wp )
   ELSE
      alpha(j) = CONJG( CMPLX( alphar(i), alphai(i), KIND = wp ) )
      END IF
END DO

IF( nA /= ZERO )THEN
   DO j = 1, k
      i = indx(j)
      CTEMP(1:n,j) = ( alpha(j)**2 )*CMPLX( V(1:n,i), V(1:n,i+1), KIND = wp )
   END DO

   TEMP(1:n,1:n) = CMPLX( A(1:n,1:n), KIND = wp )
   CALL ZGEMM( Conj, 'NoConj', n, k, n, CONE, TEMP, n, CTEMP, n, CZERO, R, n )

ELSE
   R(1:n,1:k) = CZERO
END IF

IF( nB /= ZERO )THEN
   DO j = 1, k
      i = indx(j)
      CTEMP(1:n,j) = ( alpha(j)*beta(i) )*CMPLX( V(1:n,i), V(1:n,i+1),         &
                     KIND = wp )
   END DO

   TEMP(1:n,1:n) = CMPLX( B(1:n,1:n), KIND = wp )
   CALL ZGEMM( Conj, 'NoConj', n, k, n, CONE, TEMP, n, CTEMP, n, CONE, R, n )

END IF

IF( nC /= ZERO )THEN
   DO j = 1, k
      i = indx(j)
      CTEMP(1:n,j) = ( beta(i)**2 )*CMPLX( V(1:n,i), V(1:n,i+1), KIND = wp )
   END DO

   TEMP(1:n,1:n) = CMPLX( C(1:n,1:n), KIND = wp )
   CALL ZGEMM( Conj, 'NoConj', n, k, n, CONE, TEMP, n, CTEMP, n, CONE, R, n )

END IF
DEALLOCATE( alpha, CTEMP, TEMP, STAT = iwarn_STAT )

RETURN

END SUBROUTINE DLAG3C
