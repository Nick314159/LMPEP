SUBROUTINE DLAG3V( n, m, kr, indxr, kc, indxc, alphar, alphai, beta,           &
                   A, ldA, B, ldB, C, ldC, nA, nB, nC, V1, ldV1, V2, ldV2,     &
                   V, ldV, INFO )

! Purpose
! =======

! To find the vectors with smallest backward error of a quadratic eigenvalue
! problem from two sets of n element vectors.  The vectors are the columns of
! the two n by m matrices V1 and V2 and the jth eigenvector of each is compared,
! with the eigenvector of smallest backward error being stored in V.
!
! The eigenvectors in V1 and V2 are assumed to be stored in the same way as
! required by DG3EVX.  That is, if the jth eigenvalue is real then the jth
! columns of V1 and V2 are potential jth eigenvectors.  If the jth and (j+1)th
! eigenvalues form a complex conjugate pair, then (V(:,j) + i*V(:,j+1)) and
! (V(:,j) - i*V(:,j+1)) are the corresponding complex conjugate eigenvectors.
!
! The quadratic eigenvalue problem is given in the form
!
!   ( (alpha(j)**2)*A + (alpha(j)*beta(j))*B + (beta(j)**2)*C )*x(j) = 0,
!
! where A, B and C are n by n real matrices, alpha and beta are vectors,
! with alpha = alphar + i*alphai and beta real, and x(j) is a vector.  

! Arguments
! =========

! n      - (input) INTEGER
!          On entry: the order of the matrices A, B and C.
!          Constraint: n >= 0.

! m      - (input) INTEGER
!          On entry: the number of columns of the matrices V1, V2 and V.
!          Constraint: m >= 0.

! kr    - (input) INTEGER
!         On entry: the number of columns of V1 (and V2) that represent the
!                   real vectors.
!         Constraint: kr >= 0.

! indxr - (input) INTEGER array, dimension (kr)
!         On entry: the indices of the real eigenvalues.

! kc    - (input) INTEGER
!         On entry: the number of columns of V1 (and V2) that represent the
!                   complex vectors.
!         Constraint: kc >= 0.

! indxc - (input) INTEGER array, dimension (kc)
!         On entry: the indices of the complex eigenvalues.

! alphar - (input) REAL(KIND=wp) array, dimension (m)
!          On entry: alphar(j), j = 1, 2, ..., m, contains the real parts of the
!                    elements of the vector alpha.

! alphai - (input) REAL(KIND=wp) array, dimension (m)
!          On entry: alphai(j), j = 1, 2, ..., m, contains the imaginary parts
!                    elements of the vector alpha.

! beta   - (input) REAL(KIND=wp) array, dimension (m)
!          On entry: beta(j), j = 1, 2, ..., m, contains the elements of the
!                    vector beta.

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

! V1     - (input/output) REAL(KIND=wp) array, dimension (ldV1,*)
!          Note: the second dimension of the array V1 must be at least m.
!          On entry: the n by m matrix V1.
!          On exit: each eigenvector in V1 is normalized in the 2-norm.

! ldV1   - (input) INTEGER
!          On entry: the leading dimension of the array V1.
!          Constraint: ldV1 >= n.

! V2     - (input/output) REAL(KIND=wp) array, dimension (ldV2,*)
!          Note: the second dimension of the array V2 must be at least m.
!          On entry: the n by m matrix V2.
!          On exit: each eigenvector in V2 is normalized in the 2-norm.

! ldV2   - (input) INTEGER
!          On entry: the leading dimension of the array V2.
!          Constraint: ldV2 >= n.

! V      - (output) REAL(KIND=wp) array, dimension (ldV,*)
!          Note: the second dimension of the array V must be at least m.
!          On exit: the n by m matrix V.  Each eigenvector in V is normalized in
!                   the 2-norm and the first element of largest absolute value
!                   of each vector is real and positive

! ldV    - (input) INTEGER
!          On entry: the leading dimension of the array V.
!          Constraint: ldV >= n.

! INFO    - (output) INTEGER 
!           On exit: INFO=0 unless the routine detects an error.
!           INFO > 0
!              Allocation of memory failed.  INFO returns the value of the
!              STAT  flag from the Fortran ALLOCATE statement of the compiler
!              with which the routine was compiled.

!      .. Implicit None Statement ..
       IMPLICIT NONE
!      .. Parameters ..
       INTEGER, PARAMETER             :: wp = KIND(0.0D0)
       REAL (KIND=wp), PARAMETER      :: ZERO = 0.0E+0_wp
!      .. Scalar Arguments ..
       INTEGER, INTENT (IN)           :: kc, kr, ldA, ldB, ldC, ldV, ldV1,     &
                                         ldV2, n, m
       INTEGER, INTENT (OUT)          :: INFO
       REAL (KIND=wp), INTENT (IN)    :: nA, nB, nC
!      .. Array Arguments ..
       INTEGER, INTENT (IN)           :: indxc(kc), indxr(kr)
       REAL (KIND=wp), INTENT (IN)    :: alphai(m), alphar(m), beta(m)
       REAL (KIND=wp), INTENT (INOUT) :: A(ldA,*), B(ldB, *), C(ldC,*),        &
                                         V1(ldV1,*), V2(ldV2,*)
       REAL (KIND=wp), INTENT (OUT)   :: V(ldV, *)
!      .. Local Scalars ..
       INTEGER                        :: iwarn_STAT, j, k, mk
!      .. Local Arrays ..
       COMPLEX (KIND=wp), ALLOCATABLE :: Rc(:,:)
       REAL (KIND=wp), ALLOCATABLE    :: abs1(:), abs2(:), res1(:), res2(:),   &
                                         Rr(:,:)
!      .. External Subroutines ..
       EXTERNAL                          DLAG3C, DLAG3R, DLASGE
!      .. Intrinsic Functions ..
       INTRINSIC                         ABS, CMPLX, SUM
!      .. Executable Statements ..
       CONTINUE

! Test input arguments
INFO = 0

IF( (n == 0).OR.(m == 0) )THEN
   RETURN
END IF

mk = MAX( kr, kc )
ALLOCATE( abs1(mk), abs2(mk), res1(mk), res2(mk), Rr(n,kr), STAT = iwarn_STAT )
IF( iwarn_STAT /= 0 )THEN
   INFO = iwarn_STAT
   RETURN
END IF

! Need to find the vectors with smallest backward error

! Normalize the vectors so that norm2(y) = 1
CALL DLASGE( 0, n, m, alphai, V1, ldV1, INFO )
IF( iwarn_STAT /= 0 )THEN
   INFO = iwarn_STAT
   RETURN
END IF
CALL DLASGE( 0, n, m, alphai, V2, ldV2, INFO )
IF( iwarn_STAT /= 0 )THEN
   INFO = iwarn_STAT
   RETURN
END IF

! First deal with the real case
DO j = 1, kr
   k = indxr(j)
   abs1(j) = SUM( ABS( V1(1:n,k) ) )
   abs2(j) = SUM( ABS( V2(1:n,k) ) )
END DO

! Compute the residuals corresponding to the real eigenvalues for V1
CALL DLAG3R( 'NoTrans', n, m, kr, indxr, alphar, beta, A, ldA, B, ldB, C, ldC, &
             nA, nB, nC, V1, ldV1, Rr, n, INFO )
IF( INFO > 0 )THEN
   RETURN
END IF

DO j = 1, kr
   IF( abs1(j) /= ZERO )THEN
      res1(j) = SUM( ABS( Rr(:,j) ) )/abs1(j)
   ELSE
      res1(j) = ZERO
   END IF
END DO

! Compute the residuals corresponding to the real eigenvalues for V2
CALL DLAG3R( 'NoTrans', n, m, kr, indxr, alphar, beta, A, ldA, B, ldB, C, ldC, &
             nA, nB, nC, V2, ldV2, Rr, n, INFO )
IF( INFO > 0 )THEN
   RETURN
END IF

DO j = 1, kr
   IF( abs2(j) /= ZERO )THEN
      res2(j) = SUM( ABS( Rr(:,j) ) )/abs2(j)
   ELSE
      res2(j) = ZERO
   END IF
END DO

DO j = 1, kr
   k = indxr(j)
   IF( ( abs1(j) /= ZERO ).AND.( abs2(j) /= ZERO ) )THEN
      IF( res1(j) <= res2(j) )THEN
         V(1:n,k) = V1(1:n,k)
      ELSE
         V(1:n,k) = V2(1:n,k)
      END IF
   ELSE IF( abs1(j) /= ZERO )THEN
      V(1:n,k) = V1(1:n,k)
   ELSE
      V(1:n,k) = V2(1:n,k)
   END IF
END DO

DEALLOCATE( Rr, STAT = iwarn_STAT )

! Now deal with the complex case
DO j = 1, kc
   k = indxc(j)
   abs1(j) = SUM( ABS( CMPLX( V1(1:n,k), V1(1:n,k+1), KIND = wp ) ) )
   abs2(j) = SUM( ABS( CMPLX( V2(1:n,k), V2(1:n,k+1), KIND = wp ) ) )
END DO

ALLOCATE( Rc(n,kc), STAT = iwarn_STAT )
IF( iwarn_STAT /= 0 )THEN
   INFO = iwarn_STAT
   RETURN
END IF

! Compute the residuals corresponding to the complex eigenvalues for V1
CALL DLAG3C( 'NoTrans', n, m, kc, indxc, alphar, alphai, beta,                 &
             A, ldA, B, ldB, C, ldC, nA, nB, nC, V1, ldV1, Rc, n, INFO )
IF( INFO > 0 )THEN
   RETURN
END IF

DO j = 1, kc
   IF( abs1(j) /= ZERO )THEN
      res1(j) = SUM( ABS( Rc(:,j) ) )/abs1(j)
   ELSE
      res1(j) = ZERO
   END IF
END DO

! Compute the residuals corresponding to the complex eigenvalues for V2
CALL DLAG3C( 'NoTrans', n, m, kc, indxc, alphar, alphai, beta,                 &
             A, ldA, B, ldB, C, ldC, nA, nB, nC, V2, ldV2, Rc, n, INFO )
IF( INFO > 0 )THEN
   RETURN
END IF

DO j = 1, kc
   IF( abs2(j) /= ZERO )THEN
      res2(j) = SUM( ABS( Rc(:,j) ) )/abs2(j)
   ELSE
      res2(j) = ZERO
   END IF
END DO

DO j = 1, kc
   k = indxc(j)
   IF( ( abs1(j) /= ZERO ).AND.( abs2(j) /= ZERO ) )THEN
      IF( res1(j) <= res2(j) )THEN
         V(1:n,k) = V1(1:n,k)
         V(1:n,k+1) = V1(1:n,k+1)
      ELSE
         V(1:n,k) = V2(1:n,k)
         V(1:n,k+1) = V2(1:n,k+1)
      END IF
   ELSE IF( abs1(j) /= ZERO )THEN
      V(1:n,k) = V1(1:n,k)
      V(1:n,k+1) = V1(1:n,k+1)
   ELSE
      V(1:n,k) = V2(1:n,k)
      V(1:n,k+1) = V2(1:n,k+1)
   END IF
END DO

DEALLOCATE( abs1, abs2, res1, res2, Rc, STAT = iwarn_STAT )

! Make the first element of largest absolute value of each vector real and
! positive
CALL DLASGE( 1, n, m, alphai, V, ldV, INFO )
IF( iwarn_STAT /= 0 )THEN
   INFO = iwarn_STAT
   RETURN
END IF

RETURN

END SUBROUTINE DLAG3V
