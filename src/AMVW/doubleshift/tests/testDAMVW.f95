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
  INTEGER, PARAMETER :: dp=kind(0.0D0), itmax=60

  ! input variables
  integer :: N, FLAG, NEWTNUM, rsize, maxDegree
  complex(kind(1d0)), allocatable :: COEFFS(:), ALLROOTS(:,:), ROOTS(:), WPOLY(:)
  double precision, allocatable :: POLY(:), REIGS(:), IEIGS(:), RESIDUALS(:,:)
  integer, allocatable :: ITS(:),  seed(:)
  
  ! compute variables
  integer :: ii, noits, mri, mri1, mri2, mri3, kk, i, m, j
  integer :: clock_start, clock_end, clock_rate 
  real :: time
  double precision :: rpart, ipart, temp, mr, mr1, mr2, mr3
  character(len=32) :: arg
  REAL(dp), DIMENSION(:), ALLOCATABLE :: timeStats, radStats
<<<<<<< HEAD
  
  CHARACTER(*), PARAMETER :: resultsDir="/home/thomas/Documents/FORTRAN/Nick/LMPEPtests/results/"
  !CHARACTER(*), PARAMETER :: resultsDir="/home/nsteckley/Documents/Personal/Cameron/LMPEP/tests/results/"
=======
  COMPLEX(dp) :: a, b, t
  CHARACTER(*), PARAMETER :: resultsDir="/home/thomas/Documents/FORTRAN/Nick/LMPEPtests/results/"
  !CHARACTER(*), PARAMETER :: resultsDir="/home/nsteckley/Documents/Personal/Cameron/LMPEP/tests/results/"
  REAL(dp), DIMENSION(:), ALLOCATABLE :: radius
>>>>>>> 88cb91edf1d05d0410007c86ed61fbe5475ccf00
  
  FLAG = 1
  if (iargc()>0) then
     if (iargc()>2) then
        FLAG=101
     end if
     call RANDOM_SEED(size = rsize)
     allocate(seed(rsize))
     call RANDOM_SEED(GET = seed)

     !print*, iargc()
     !print*, arg
     if (iargc()>1) then
        call getarg(2, arg)
        print*, arg

        read (arg,'(I10)') ii
        !print*, seed
        seed(1) = ii
        seed(2) = ii+1000000
        seed(3) = ii+2100000
        seed(4) = ii+3210000
        seed(5) = ii+43210000
        seed(6) = ii+5432100
        seed(7) = ii+6543210
        seed(8) = ii+7654321
        seed(9) = ii+8765432
        seed(10) = ii+9876543
        seed(11) = ii+10987654
        seed(12) = ii+11109876
        call RANDOM_SEED(PUT = seed)
     else
        call init_random_seed()
     end if

     call getarg(1, arg)
     read (arg,'(I10)') kk
  else
     kk = 4096
  end if
 
  maxDegree = kk
  N = 10
  m = 10
  ALLOCATE(timeStats(m), radStats(m))

  NEWTNUM = 1
  
   OPEN(UNIT=1,FILE=resultsDir//"outputAMVW.csv")
    WRITE(1, '(A)',  advance='no') 'DEGREE,       '
    WRITE(1, '(A)',  advance='no') 'AMVW TIME,    '
    WRITE(1, '(A)',  advance='no') 'AMVW Radius,  '
    WRITE(1, *)
    
    
  DO WHILE (N < maxDegree)
    WRITE(1, '(i6)', advance='no') N
    WRITE(1, '(A)', advance='no') ', '
    DO i=1,m
  
      allocate(POLY(N),REIGS(N),IEIGS(N),ITS(N),COEFFS(1),ALLROOTS(N,NEWTNUM+1),RESIDUALS(N,3*(NEWTNUM+1)))
      allocate(WPOLY(N),ROOTS(N))
  
      ALLOCATE(radius(1:N))
      call DNORMALPOLY(N,POLY)
 
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
      !radStats(i)=MAXVAL(radius)
      deallocate(POLY,REIGS,IEIGS,ITS,COEFFS,ALLROOTS,RESIDUALS,WPOLY,ROOTS);
      deallocate(radius)
    end do 
    WRITE(1,'(20G15.4)', advance='no') SUM(timeStats)/10
    WRITE(1, '(A)', advance='no') ', '
    WRITE(1,'(20G15.4)', advance='no') SUM(radStats)/10
    WRITE(1, *)
    N = 2 * N
  ENDDO
  CLOSE(UNIT=1)
END PROGRAM

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
  a=p(d+1)
  DO k=d,1,-1
    a=t*a+p(k)
  ENDDO
ELSEIF(der==1) THEN
  a=d*p(d+1)
  DO k=d,2,-1
    a=t*a+(k-1)*p(k)
  ENDDO
ELSE
  a=d*(d-1)*p(d+1)
  DO k=d,3,-1
    a=t*a+(k-1)*(k-2)*p(k)
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
END SUBROUTINE zseval
