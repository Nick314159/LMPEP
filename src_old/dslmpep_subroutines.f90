MODULE dslmpep_subroutines

IMPLICIT NONE
!parameters
INTEGER, PARAMETER :: dp=KIND(0.0D0), itmax=50
REAL(dp), PARAMETER :: eps=EPSILON(0.0_dp), big=HUGE(0.0_dp)
REAL(dp), PARAMETER :: zero=0.0_dp, one=1.0_dp
COMPLEX(dp), PARAMETER :: czero=DCMPLX(zero), cone=DCMPLX(one)

CONTAINS

!************************************************************************
!			SUBROUTINE DSLCORR				*
!************************************************************************
! Compute the Laguerre correction term for scalar polynomial p with real*
! coeffs of degree d at real number (tr,0). The magnitude of the	*
! coeffs are stored in alpha, all n previously found roots are stored	*
! in (er,ei). Result is returned in complex number (tr,ti), and backward*
! error of current approximation is stored in be.			*
!************************************************************************
SUBROUTINE dslcorr(p, alpha, er, ei, berr, tol, check, d, i)
IMPLICIT NONE
!scalar arguments
LOGICAL, INTENT(INOUT) :: check
INTEGER, INTENT(IN) :: d, i
REAL(dp), INTENT(IN) :: tol
REAL(dp), INTENT(INOUT) :: berr
!array arguments
REAL(dp), INTENT(IN) :: p(*), alpha(*)
REAL(dp), INTENT(INOUT) :: er(*), ei(*)
!local scalars
INTEGER :: k
REAL(dp) :: a, b, c, t
DOUBLE COMPLEX :: x1, x2, y1, y2
!intrinsic procedures
INTRINSIC :: DABS, DBLE, DCMPLX, DIMAG, ZABS, ZSQRT

!initiate variables
x1=czero; x2=czero
DO k=1,i-1
  y1=1/DCMPLX(er(i)-er(k),-ei(k))
  x1=x1+y1
  x2=x2+y1**2
ENDDO
t=er(i)
!split into 2 cases
IF(DABS(t)>1) THEN
  !compute a=revp, berr
  t=1/t
  CALL drevseval(p, t, a, d, 0)
  CALL drevseval(alpha, DABS(t), berr, d, 0)
  berr=MIN(DABS(a)/berr, DABS(a))
  IF(berr<eps) THEN
    ei(i)=zero
    check=.TRUE.
    RETURN
  ELSE
    !compute b=revp', c=revp''
    CALL drevseval(p, t, b, d, 1)
    CALL drevseval(p, t, c, d, 2)
    !compute y1=p'/p and y2=(p'/p)'
    y1=t*(d-t*(b/a))
    y2=t**2*(d-2*t*(b/a)+t**2*((b/a)**2-c/a))
  ENDIF
ELSE
  !compute a=p, berr
  CALL dseval(p, t, a, d, 0)
  CALL dseval(alpha, DABS(t), berr, d, 0)
  berr=MIN(DABS(a)/berr, DABS(a))
  IF(berr<eps) THEN
    ei(i)=zero
    check=.TRUE.
    RETURN
  ELSE
    !compute b=p', c=p''
    CALL dseval(p, t, b, d, 1)
    CALL dseval(p, t, c, d, 2)
    !compute y1=p'/p and y2=(p'/p)'
    y1=b/a
    y2=y1**2-c/a
  ENDIF
ENDIF
!remove previously found roots
x1=y1-x1
x2=y2-x2
k=d-i+1
!denominator of Laguerre correction term
y1=ZSQRT((k-1)*(k*x2-x1**2))
y2=x1-y1; y1=x1+y1
IF(ZABS(y1)>=ZABS(y2)) THEN
  y1=k/y1
  IF(ZABS(y1)<tol) THEN
    ei(i)=zero
    check=.TRUE.
  ELSE
    er(i)=er(i)-DBLE(y1)
    ei(i)=-DIMAG(y1)
  ENDIF
ELSE
  y2=k/y2
  IF(ZABS(y2)<tol) THEN
    ei(i)=zero
    check=.TRUE.
  ELSE
    er(i)=er(i)-DBLE(y2)
    ei(i)=-DIMAG(y2)
  ENDIF
ENDIF
RETURN
END SUBROUTINE dslcorr

!************************************************************************
!			SUBROUTINE ZSLCORR				*
!************************************************************************
! Compute the Laguerre correction term for scalar polynomial p with real*
! coeffs of degree d at complex number (tr,ti). The magnitude of the	*
! coeffs are stored in alpha, all n previously found roots are stored	*
! in (er,ei). Result is return in complex number (tr,ti), and backward	*
! error of current approximation is stored in be.			*
!************************************************************************
SUBROUTINE zslcorr(p, alpha, er, ei, berr, tol, check, d, i)
IMPLICIT NONE
!scalar arguments
LOGICAL, INTENT(INOUT) :: check
INTEGER, INTENT(IN) :: d, i
REAL(dp), INTENT(IN) :: tol
REAL(dp), INTENT(INOUT) :: berr
!array arguments
REAL(dp), INTENT(IN) :: p(*), alpha(*)
REAL(dp), INTENT(INOUT) :: er(*), ei(*)
!local scalars
INTEGER :: k
COMPLEX(dp) :: a, b, c, t, x1, x2, y1, y2
!intrinsic procedures
INTRINSIC :: DBLE, DCMPLX, DIMAG, ZABS, ZSQRT

!initiate variables
x1=czero; x2=czero
DO k=1,i-1
  y1=1/DCMPLX(er(i)-er(k),ei(i)-ei(k))
  x1=x1+y1
  x2=x2+y1**2
ENDDO
t=DCMPLX(er(i),ei(i))
!split into 2 cases
IF(ZABS(t)>1) THEN
  !compute a=revp, berr
  t=1/t
  CALL zrevseval(p, t, a, d, 0)
  CALL drevseval(alpha, ZABS(t), berr, d, 0)
  berr=MIN(ZABS(a)/berr, ZABS(a))
  IF(berr<eps) THEN
    check=.TRUE.
    RETURN
  ELSE
    !compute b=revp', c=revp''
    CALL zrevseval(p, t, b, d, 1)
    CALL zrevseval(p, t, c, d, 2)
    !compute y1=p'/p and y2=(p'/p)'
    y1=t*(d-t*(b/a))
    y2=t**2*(d-2*t*(b/a)+t**2*((b/a)**2-c/a))
  ENDIF
ELSE
  !compute a=p, berr
  CALL zseval(p, t, a, d, 0)
  CALL dseval(alpha, ZABS(t), berr, d, 0)
  berr=MIN(ZABS(a)/berr, ZABS(a))
  IF(berr<eps) THEN
    check=.TRUE.
    RETURN
  ELSE
    !solve for b=p', c=p''
    CALL zseval(p, t, b, d, 1)
    CALL zseval(p, t, c, d, 2)
    !compute y1=p'/p and y2=(p'/p)'
    y1=b/a
    y2=y1**2-c/a
  ENDIF
ENDIF
!remove previously found roots
x1=y1-x1
x2=y2-x2
k=d-i+1
!denominator of Laguerre correction term
y1=ZSQRT((k-1)*(k*x2-x1**2))
y2=x1-y1; y1=x1+y1
IF(ZABS(y1)>=ZABS(y2)) THEN
  y1=k/y1
  IF(ZABS(y1)<tol) THEN
    check=.TRUE.
  ELSE
    er(i)=er(i)-DBLE(y1)
    ei(i)=ei(i)-DIMAG(y1)
  ENDIF
ELSE
  y2=k/y2
  IF(ZABS(y2)<tol) THEN
    check=.TRUE.
  ELSE
    er(i)=er(i)-DBLE(y2)
    ei(i)=ei(i)-DIMAG(y2)
  ENDIF
ENDIF
RETURN
END SUBROUTINE zslcorr

!************************************************************************
!			SUBROUTINE DREVSEVAL				*
!************************************************************************
! Evaluate reversal of scalar polynomial p with real coeffs of degree d,*
! and its der=0,1,2 derivatives at real number 1/t, where |t|>1. 	*
! Returns evaluation in a.						*
!************************************************************************
SUBROUTINE drevseval(p, t, a, d, der)
IMPLICIT NONE
!scalar arguments
INTEGER, INTENT(IN) :: d, der
REAL(dp), INTENT(IN) :: t
REAL(dp), INTENT(INOUT) :: a
!array arguments
REAL(dp), INTENT(IN) :: p(*)
!local scalars
INTEGER :: k

IF(der==0) THEN
  a=p(1)
  DO k=2,d+1
    a=t*a+p(k)
  ENDDO
ELSEIF(der==1) THEN
  a=d*p(1)
  DO k=2,d
    a=t*a+(d-k+1)*p(k)
  ENDDO
ELSE
  a=d*(d-1)*p(1)
  DO k=2,d-1
    a=t*a+(d-k+1)*(d-k)*p(k)
  ENDDO
ENDIF
RETURN
END SUBROUTINE drevseval



!************************************************************************
!			SUBROUTINE ZREVSEVAL				*
!************************************************************************
! Evaluate reversal of scalar polynomial p with real coeffs of degree d,*
! and its der=0,1,2 derivatives at complex number 1/t, where |t|>1. 	*
! Returns evaluation in a.						*
!************************************************************************
SUBROUTINE zrevseval(p, t, a, d, der)
IMPLICIT NONE
!scalar arguments
INTEGER, INTENT(IN) :: d, der
COMPLEX(dp), INTENT(IN) :: t
COMPLEX(dp), INTENT(INOUT) :: a
!array arguments
REAL(dp), INTENT(IN) :: p(*)
!local scalars
INTEGER :: k

IF(der==0) THEN
  a=p(1)
  DO k=2,d+1
    a=t*a+p(k)
  ENDDO
ELSEIF(der==1) THEN
  a=d*p(1)
  DO k=2,d
    a=t*a+(d-k+1)*p(k)
  ENDDO
ELSE
  a=d*(d-1)*P(1)
  DO k=2,d-1
    a=t*a+(d-k+1)*(d-k)*p(k)
  ENDDO
ENDIF
RETURN
END SUBROUTINE zrevseval


!************************************************************************
!			SUBROUTINE DSSTART				*
!************************************************************************
! Compute the initial esimates of the roots of a real scalar polynomial *
! using the Newton Polygon method. The absolute value of the 		*
! coefficients of the polynomial of degree d are stored in alpha, the 	*
! real and imaginary parts of the initial estiamtes are returned in er 	*
! and ei.         							*
!************************************************************************
SUBROUTINE dsstart(alpha, er, ei, d)
IMPLICIT NONE
!scalar arguments
INTEGER, INTENT(IN) :: d
!array arguments
REAL(dp), INTENT(IN) :: alpha(*)
REAL(dp), INTENT(INOUT) :: er(*), ei(*)
!parameters
REAL(dp), PARAMETER :: pi2 = 6.2831853071795865_dp, sigma = 0.7_dp
!local scalars
INTEGER :: c, i, iold, j, nzeros
REAL(dp) :: ang, r, th
!local arrays
LOGICAL, DIMENSION(d+1) :: h
REAL(dp), DIMENSION(d+1) :: a
!intrinsic procedures
INTRINSIC :: DCOS, DEXP, DLOG, DSIN

!compute log(alpha)
DO i=1,d+1
  IF(alpha(i)>=eps) THEN
    a(i)=DLOG(alpha(i))
  ELSE
    a(i)=-one
  ENDIF
ENDDO
!compute upper convex hull
CALL cnvex(d+1,a,h)
!compute initial estimates
iold=1; c=0; th=pi2/d
DO i=2,d+1
  IF(h(i)) THEN
    nzeros=i-iold
    r=DEXP((a(iold)-a(i))/nzeros)
    ang=pi2/nzeros
    DO j=1,nzeros
      er(c+j)=r*DCOS(ang*j+th*i+sigma)
      ei(c+j)=r*DSIN(ang*j+th*i+sigma)
    ENDDO
    c=c+nzeros
    iold=i
  ENDIF
ENDDO
RETURN
END SUBROUTINE dsstart

!************************************************************************
!			SUBROUTINE DCMOD				*
!************************************************************************
! Compute the module of the complex number a+bi, while avoiding harmul  *
! overflow and underflow. Sqrt(a^(2)+b^(2))				*
!************************************************************************
FUNCTION dcmod(a, b) RESULT(r)
IMPLICIT NONE
!scalar arguments
REAL(dp), INTENT(IN) :: a, b
!local scalars
REAL(dp) :: r

IF(DABS(a)<eps .AND. DABS(b)<eps) THEN
  r=zero
  RETURN
ENDIF

IF(DABS(a)<DABS(b)) THEN
  r=DABS(b)*DSQRT(1+(a/b)**2)
ELSE
  r=DABS(a)*DSQRT(1+(b/a)**2)
ENDIF
RETURN
END FUNCTION dcmod

!************************************************************************
!                             SUBROUTINE CNVEX                          *
!************************************************************************
! Compute  the upper convex hull of the set (i,a(i)), i.e., the set of  *
! vertices (i_k,a(i_k)), k=1,2,...,m, such that the points (i,a(i)) lie *
! below the straight lines passing through two consecutive vertices.    *
! The abscissae of the vertices of the convex hull equal the indices of *
! the TRUE  components of the logical output vector H.                  *
! The used method requires O(nlog n) comparisons and is based on a      *
! divide-and-conquer technique. Once the upper convex hull of two       *
! contiguous sets  (say, {(1,a(1)),(2,a(2)),...,(k,a(k))} and           *
! {(k,a(k)), (k+1,a(k+1)),...,(q,a(q))}) have been computed, then       *
! the upper convex hull of their union is provided by the subroutine    *
! CMERGE. The program starts with sets made up by two consecutive       *
! points, which trivially constitute a convex hull, then obtains sets   *
! of 3,5,9... points,  up to  arrive at the entire set.                 *
! The program uses the subroutine  CMERGE; the subroutine CMERGE uses   *
! the subroutines LEFT, RIGHT and CTEST. The latter tests the convexity *
! of the angle formed by the points (i,a(i)), (j,a(j)), (k,a(k)) in the *
! vertex (j,a(j)) up to within a given tolerance TOLER, where i<j<k.    *
!************************************************************************
SUBROUTINE cnvex(n, a, h)
IMPLICIT NONE
INTEGER, INTENT(IN)        :: n
LOGICAL, INTENT(OUT)       :: h(:)
REAL(dp), INTENT(IN) :: a(:)

! Local variables
INTEGER :: i, j, k, m, nj, jc

h(1:n) = .true.

! compute K such that N-2 <= 2**K < N-1
k = INT(LOG(n-2.0_dp)/LOG(2.0_dp))
IF(2**(k+1) <= (n-2)) k = k+1

! For each M=1,2,4,8,...,2**K, consider the NJ pairs of consecutive
! sets made up by M+1 points having the common vertex
! (JC,A(JC)), where JC=M*(2*J+1)+1 and J=0,...,NJ,
! NJ = MAX(0, INT((N-2-M)/(M+M))).
! Compute the upper convex hull of their union by means of subroutine CMERGE
m = 1
DO i = 0, k
  nj = MAX(0, INT((n-2-m)/(m+m)))
  DO j = 0, nj
    jc = (j+j+1)*m+1
    CALL cmerge(n, a, jc, m, h)
  ENDDO
  m = m+m
ENDDO
RETURN
END SUBROUTINE cnvex

!************************************************************************
!                             SUBROUTINE LEFT                           *
!************************************************************************
! Given as input the integer I and the vector H of logical, compute the *
! the maximum integer IL such that IL<I and H(IL) is TRUE.              *
!************************************************************************
! Input variables:                                                      *
!     H   : vector of logical                                           *
!     I   : integer                                                     *
!************************************************************************
! Output variable:                                                      *
!     IL  : maximum integer such that IL<I, H(IL)=.TRUE.                *
!************************************************************************
SUBROUTINE left(h, i, il)
IMPLICIT NONE
INTEGER, INTENT(IN)  :: i
INTEGER, INTENT(OUT) :: il
LOGICAL, INTENT(IN)  :: h(:)

DO il = i-1, 0, -1
  IF (h(il)) RETURN
ENDDO
RETURN
END SUBROUTINE left

!************************************************************************
!                             SUBROUTINE RIGHT                          *
!************************************************************************
!************************************************************************
! Given as input the integer I and the vector H of logical, compute the *
! the minimum integer IR such that IR>I and H(IL) is TRUE.              *
!************************************************************************
!************************************************************************
! Input variables:                                                      *
!     N   : length of the vector H                                      *
!     H   : vector of logical                                           *
!     I   : integer                                                     *
!************************************************************************
! Output variable:                                                      *
!     IR  : minimum integer such that IR>I, H(IR)=.TRUE.                *
!************************************************************************
SUBROUTINE right(n, h, i, ir)
IMPLICIT NONE
INTEGER, INTENT(IN)  :: n, i
INTEGER, INTENT(OUT) :: ir
LOGICAL, INTENT(IN)  :: h(:)

DO ir = i+1, n
  IF (h(ir)) RETURN
ENDDO
RETURN
END SUBROUTINE right

!************************************************************************
!                             SUBROUTINE CMERGE                         *
!************************************************************************
! Given the upper convex hulls of two consecutive sets of pairs         *
! (j,A(j)), compute the upper convex hull of their union                *
!************************************************************************
! Input variables:                                                      *
!     N    : length of the vector A                                     *
!     A    : vector defining the points (j,A(j))                        *
!     I    : abscissa of the common vertex of the two sets              *
!     M    : the number of elements of each set is M+1                  *
!************************************************************************
! Input/Output variable:                                                *
!     H    : vector defining the vertices of the convex hull, i.e.,     *
!            H(j) is .TRUE. if (j,A(j)) is a vertex of the convex hull  *
!            This vector is used also as output.                        *
!************************************************************************
SUBROUTINE cmerge(n, a, i, m, h)
IMPLICIT NONE
INTEGER, INTENT(IN) :: n, m, i
LOGICAL, INTENT(IN OUT) :: h(:)
REAL(dp), INTENT(IN) :: a(:)

! Local variables
INTEGER :: ir, il, irr, ill
LOGICAL :: tstl, tstr

! at the left and the right of the common vertex (I,A(I)) determine
! the abscissae IL,IR, of the closest vertices of the upper convex
! hull of the left and right sets, respectively
CALL left(h, i, il)
CALL right(n, h, i, ir)

! check the convexity of the angle formed by IL,I,IR
IF (ctest(a, il, i, ir)) THEN
  RETURN
ELSE
! continue the search of a pair of vertices in the left and right
! sets which yield the upper convex hull
  h(i) = .false.
  DO
    IF (il == (i-m)) THEN
      tstl = .true.
    ELSE
      CALL left(h, il, ill)
      tstl = ctest(a, ill, il, ir)
    ENDIF
    IF (ir == MIN(n, i+m)) THEN
      tstr = .true.
    ELSE
      CALL right(n, h, ir, irr)
      tstr = ctest(a, il, ir, irr)
    ENDIF
    h(il) = tstl
    h(ir) = tstr
    IF (tstl.AND.tstr) RETURN
    IF(.NOT.tstl) il = ill
    IF(.NOT.tstr) ir = irr
  ENDDO
ENDIF
RETURN
END SUBROUTINE cmerge

!************************************************************************
!                             FUNCTION CTEST                            *
!************************************************************************
! Test the convexity of the angle formed by (IL,A(IL)), (I,A(I)),       *
! (IR,A(IR)) at the vertex (I,A(I)), up to within the tolerance         *
! TOLER. If convexity holds then the function is set to .TRUE.,         *
! otherwise CTEST=.FALSE. The parameter TOLER is set to 0.4 by default. *
!************************************************************************
! Input variables:                                                      *
!     A       : vector of double                                        *
!     IL,I,IR : integers such that IL < I < IR                          *
!************************************************************************
! Output:                                                               *
!     .TRUE. if the angle formed by (IL,A(IL)), (I,A(I)), (IR,A(IR)) at *
!            the vertex (I,A(I)), is convex up to within the tolerance  *
!            TOLER, i.e., if                                            *
!            (A(I)-A(IL))*(IR-I)-(A(IR)-A(I))*(I-IL)>TOLER.             *
!     .FALSE.,  otherwise.                                              *
!************************************************************************
FUNCTION ctest(a, il, i, ir) RESULT(OK)
IMPLICIT NONE
INTEGER, INTENT(IN) :: i, il, ir
REAL(dp), INTENT(IN) :: a(:)
LOGICAL :: OK

! Local variables
REAL(dp) :: s1, s2
REAL(dp), PARAMETER :: tol = 0.4_dp

s1 = a(i) - a(il)
s2 = a(ir) - a(i)
s1 = s1*(ir-i)
s2 = s2*(i-il)
OK = .false.
IF(s1 > (s2+tol)) OK = .true.
RETURN
END FUNCTION ctest

END MODULE dslmpep_subroutines
