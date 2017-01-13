!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Aurentz² Mach³ Vandebril³ Watkins²
!
! ²Dept. Mathematics, Washington State University
! ³Dept. Computer Science, KU Leuven
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Last modified 22 August 2014
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! 
! Test program for random polynomials of given 
! size
!
! The roots are computed by Francis's implicitly
! shifted QR algorithm via the Companion matrix.
! The rank-structure in the iterates is used.
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! input parameter
!
! 1) problem size, default 4096
!
! 2) seed for random number generator, 
!    default random seed based on CPU clock, 
!    set fixed seed for reproducibility
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
program testDAMVW

  implicit none
  INTEGER, PARAMETER :: dp=KIND(0.0D0)

  ! input variables
  integer :: N, FLAG, NEWTNUM, rsize, maxDegree, jumpSize
  complex(kind(1d0)), allocatable :: COEFFS(:), ALLROOTS(:,:), ROOTS(:), WPOLY(:)
  double precision, allocatable :: POLY(:), REIGS(:), IEIGS(:), RESIDUALS(:,:)
  integer, allocatable :: ITS(:),  seed(:)
  ! compute variables
  integer :: i, m, j
  integer :: clock_start, clock_end, clock_rate 
  character(len=32) :: arg
  REAL(dp), DIMENSION(:), ALLOCATABLE :: timeStats, radStats
  
CHARACTER(*), PARAMETER :: resultsDir="/home/nsteckley/Documents/Personal/Cameron/LMPEP/results/"
  COMPLEX(dp) :: a, b, t
  REAL(dp), DIMENSION(:), ALLOCATABLE :: radius

  
  FLAG = 1

  CALL RANDOM_SEED(size = rsize)
  ALLOCATE(seed(rsize))
  CALL RANDOM_SEED(GET = seed)
  CALL init_random_seed()
  CALL GETARG(1, arg)
  READ (arg,'(I10)') N
  CALL GETARG(2, arg)
  READ (arg,'(I10)') maxDegree
  CALL GETARG(3, arg)
  READ (arg,'(I10)') m
  !CALL GETARG(3, arg)
  !READ (arg,'(I10)') jumpSize 
  jumpSize = 2

  ALLOCATE(timeStats(m), radStats(m))

  NEWTNUM = 1
  
  OPEN(UNIT=1,FILE=resultsDir//"outputAMVW.csv")
  WRITE(1, '(A)',  advance='no') 'DEGREE,      '
  WRITE(1, '(A)',  advance='no') 'AMVW TIME,       '
  WRITE(1, '(A)',  advance='no') 'AMVW Radius  '
  WRITE(1, *)
      
  DO WHILE (N < maxDegree)
    WRITE(1, '(i6)', advance='no') N
    WRITE(1, '(A)', advance='no') ', '
    DO i=1,m
  
      ALLOCATE(POLY(N),REIGS(N),IEIGS(N),ITS(N),COEFFS(1),ALLROOTS(N,NEWTNUM+1),RESIDUALS(N,3*(NEWTNUM+1)))
      ALLOCATE(WPOLY(N),ROOTS(N))
  
      ALLOCATE(radius(1:N))
      CALL DNORMALPOLY(N,POLY)
      ! start timer
      CALL SYSTEM_CLOCK(COUNT_RATE=clock_rate) 
      CALL SYSTEM_CLOCK(COUNT=clock_start)
      ! compute roots
      call DAMVW(N,POLY,REIGS,IEIGS,ITS,FLAG)
      ! stop timer
      CALL SYSTEM_CLOCK(COUNT=clock_end)  
      timeStats(i) = DBLE(clock_end-clock_start)/DBLE(clock_rate)
  
      DO j=1,N
        t=DCMPLX(REIGS(j),IEIGS(j))
        IF(ZABS(t)>1) THEN
          t=1/t
          CALL zrevseval(POLY, t, a, N, 0)
          CALL zrevseval(POLY, t, b, N, 1)
          radius(j)=ZABS(a)/ZABS(N*a-b*t)
        ELSE
          CALL zseval(POLY, t, a, N, 0)
          CALL zseval(POLY, t, b, N, 1)
          radius(j)=ZABS(a)/(ZABS(t)*ZABS(b))
        ENDIF
      ENDDO
      radStats(i)=MAXVAL(radius)
      DEALLOCATE(POLY,REIGS,IEIGS,ITS,COEFFS,ALLROOTS,RESIDUALS,WPOLY,ROOTS);
      DEALLOCATE(radius)
    END DO 
    WRITE(1,'(20G15.4)', advance='no') SUM(timeStats)/m
    WRITE(1, '(A)', advance='no') ', '
    WRITE(1,'(20G15.4)', advance='no') SUM(radStats)/m
    WRITE(1, *)
    N = jumpSize * N
  ENDDO
  CLOSE(UNIT=1)
CONTAINS

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
  a=p(d)
  DO k=d-1,1,-1
    a=t*a+p(k)
  ENDDO
  a=t*a+1
ELSEIF(der==1) THEN
  a=d*p(d)
  DO k=d-1,1,-1
    a=t*a+k*p(k)
  ENDDO
ELSE
  a=d*(d-1)*p(d)
  DO k=d-1,2,-1
    a=t*a+k*(k-1)*p(k)
  ENDDO
ENDIF
RETURN
END SUBROUTINE zrevseval

!************************************************************************
!			SUBROUTINE ZSEVAL				*
!************************************************************************
! Evaluate scalar polynomial p with real coeffs of degree d, and its	*
! der=0,1,2 derivatives at complex number t, where |t|<=1. Returns 	*
! evaluation in a.							*
!************************************************************************
SUBROUTINE zseval(p, t, a, d, der)
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
  a=1
  DO k=1,d
    a=t*a+p(k)
  ENDDO
ELSEIF(der==1) THEN
  a=d
  DO k=1,d-1
    a=t*a+(d-k)*p(k)
  ENDDO
ELSE
  a=d*(d-1)
  DO k=1,d-2
    a=t*a+(d-k)*(d-k-1)*p(k)
  ENDDO
ENDIF
RETURN
END SUBROUTINE zseval

END PROGRAM testDAMVW
