! Demo program for the computation of the eigenvalues of the 
! tridiagonal quadratic eigenvalue problems using the Ehrlich-Aberth
! and the Laguerre method
! 
! Reference: B. Plestenjak, Numerical methods for the tridiagonal
! hyperbolic quadratic eigenvalue problem, Preprint 967, IMFM, 
! Ljubljana, 2005
!
! Bor Plestenjak
! bor.plestenjak@fmf.uni-lj.si
! Department of Mathematics
! University of Ljubljana
! November 2005

PROGRAM main

USE environment
USE dgtlmpep_subroutines
USE qep3dea
USE qep3deacx
USE qep3dlag

IMPLICIT NONE

INTEGER:: mode, n, i, neg, detsgn, k, m, iseed
INTEGER, PARAMETER:: ATTEMPTS = 200
REAL(dp), ALLOCATABLE, DIMENSION(:):: a, b, c, au, bu, cu, al, bl, cl, d
REAL(dp), ALLOCATABLE, DIMENSION(:):: z
COMPLEX(dp), ALLOCATABLE, DIMENSION(:):: zcx
REAL(4):: T1
INTEGER:: mxit, iter, itermx, imax

REAL(dp), DIMENSION(:), ALLOCATABLE :: berr, cond, ferr, er, ei, ncoeff, work, x
INTEGER :: clock, clock_rate, clock_start, clock_stop
REAL(dp), DIMENSION(:,:), ALLOCATABLE :: pdl, pd, pdu, xr, xi, yr, yi
CHARACTER (LEN=64), DIMENSION(19) :: tests
REAL(dp), DIMENSION(:), ALLOCATABLE :: timeStats
ALLOCATE(timeStats(2))

tests(1) = 'data_Ex102_100_EAC.dat'
tests(2) = 'data_Ex102_100_EAR.dat'
tests(3) = 'data_Ex102_100_LAG.dat'
tests(4) = 'data_Ex102_200_EAC.dat'
tests(5) = 'data_Ex102_200_EAR.dat'
tests(6) = 'data_Ex102_200_LAG.dat'
tests(7) = 'data_Ex102_400_EAC.dat'
tests(8) = 'data_Ex102_400_EAR.dat'
tests(9) = 'data_Ex102_400_LAG.dat'
tests(10) = 'data_Ex102_800_EAC.dat'
tests(11) = 'data_Ex102_800_EAR.dat'
tests(12) = 'data_Ex102_800_LAG.dat'
tests(13) = 'data_Ex103A_100_EAC.dat'
tests(14) = 'data_Ex103A_200_EAC.dat'
tests(15) = 'data_Ex103A_400_EAC.dat'
tests(16) = 'data_Ex103B_100_EAC.dat'
tests(17) = 'data_Ex103B_200_EAC.dat'
tests(18) = 'data_Ex103B_400_EAC.dat'
tests(19) = 'data_Ex104_EAC.dat'

!Create iseed
CALL SYSTEM_CLOCK(COUNT=clock)
CALL srand(clock)
DO i=1,4
iseed(i)=MOD(irand(),4095)
ENDDO
IF(MOD(iseed(4),2)==0) THEN
iseed(4)=iseed(4)+1
ENDIF


OPEN(UNIT=1,FILE=resultsDir//"outputTri.csv")
WRITE(1, '(A)',  advance='no') 'Problem,        '
WRITE(1, '(A)',  advance='no') 'DGTLMPEP TIME,   '
WRITE(1, '(A)',  advance='no') 'QEP3D TIME'
WRITE(1, *)
d =3 
DO k=1,19
  OPEN(UNIT=2,FILE=sampleProblemsDir//tests(k))
  WRITE(1, '(A)', advance='no') tests(k)
  WRITE(1, '(A)', advance='no') ','
  READ(2,'(I2)') mode
  READ(2, '(I10)') n
  ALLOCATE(a(n), b(n), c(n), au(n), al(n), bu(n), bl(n), cu(n), cl(n))
  ALLOCATE(pdl(n-1,d+1), pd(n,d+1), pdu(n-1,d+1))
  ALLOCATE(z(2*n), zcx(2*n))
  IF (mode<5) THEN
    READ(2,*) ( a(i), i=1,n)
    READ(2,*) ( au(i), i=1,n-1)
    READ(2,*) ( b(i), i=1,n)
    READ(2,*) ( bu(i), i=1,n-1)
    READ(2,*) ( c(i), i=1,n)
    READ(2,*) ( cu(i), i=1,n-1)
 
    pd(:, 1) = a
    pdu(:, 1) = au
    pd(:, 2) = b
    pdu(:, 2) = bu
    pd(:, 3) = c
    pdu(:, 3) = cu
    
  ELSE
    READ(2,*) ( a(i), i=1,n)
    READ(2,*) ( au(i), i=1,n-1)
    READ(2,*) ( al(i), i=1,n-1)
    READ(2,*) ( b(i), i=1,n)
    READ(2,*) ( bu(i), i=1,n-1)
    READ(2,*) ( bl(i), i=1,n-1)
    READ(2,*) ( c(i), i=1,n)
    READ(2,*) ( cu(i), i=1,n-1)
    READ(2,*) ( cl(i), i=1,n-1)
    
    pd(:, 1) = a
    pdu(:, 1) = au
    pdl(:, 1) = al
    pd(:, 2) = b
    pdu(:, 2) = bu
    pdl(:, 2) = bl
    pd(:, 3) = c
    pdu(:, 3) = cu
    pdl(:, 3) = cl

  END IF
  CLOSE(UNIT=2)

  z = zero
  zcx = zero

  mxit = 400! maximal number of iteration
  iter = 0
  itermx = 500
  T1 = SECNDS(0.0)
  ! Real Ehrlich-Aberth
  IF (mode>=1 .AND. mode<=3) THEN
    CALL SYSTEM_CLOCK(count_rate=clock_rate)
    CALL SYSTEM_CLOCK(COUNT=clock_start)
    CALL reigen(a,au,b,bu,c,cu,n,z,mxit,iter,itermx,imax,mode)
    CALL SYSTEM_CLOCK(COUNT=clock_stop)
    timeStats(2) = DBLE(clock_stop-clock_start)/DBLE(clock_rate)
  END IF

  ! Laguerre
  IF (mode==4) THEN
    CALL SYSTEM_CLOCK(count_rate=clock_rate)
    CALL SYSTEM_CLOCK(COUNT=clock_start)
    CALL reigenl(a,au,b,bu,c,cu,n,z,mxit,iter,itermx,imax) 
    CALL SYSTEM_CLOCK(COUNT=clock_stop)  
    timeStats(2) = DBLE(clock_stop-clock_start)/DBLE(clock_rate)
  END IF

  ! Complex Ehrlich Aberth
  IF (mode>=5 .AND. mode<=7) THEN
    CALL SYSTEM_CLOCK(count_rate=clock_rate)
    CALL SYSTEM_CLOCK(COUNT=clock_start)
    CALL reigencx(a,au,al,b,bu,bl,c,cu,cl,n,zcx,mxit,iter,itermx,imax,mode-4)
    CALL SYSTEM_CLOCK(COUNT=clock_stop)  
    timeStats(2) = DBLE(clock_stop-clock_start)/DBLE(clock_rate)
  END IF
  
  !solve dgtlmpep
  ALLOCATE(xr(n,n*d), xi(n,n*d), yr(n,n*d), yi(n,n*d))
  ALLOCATE(berr(n*d), er(n*d), ei(n*d), ncoeff(d+1), cond(n*d), ferr(n*d))
  DO i=1,d+1
    ncoeff(i)=dlangt('F',n,pdl(1,i),pd(1,i),pdu(1,i))
  ENDDO
  CALL SYSTEM_CLOCK(count_rate=clock_rate)
  CALL SYSTEM_CLOCK(COUNT=clock_start)
  CALL dgtlm(pdl, pd, pdu, xr, xi, yr, yi, er, ei, berr, ncoeff, iseed, d, n)
  CALL SYSTEM_CLOCK(COUNT=clock_stop)
  timeStats(1) = DBLE(clock_stop-clock_start)/DBLE(clock_rate)
   
  !deallocate
  DEALLOCATE(a, b, c, au, bu, cu, al, bl, cl, z, zcx)
  DEALLOCATE(xr, xi, yr, yi)
  DEALLOCATE(berr, cond, ferr, er, ei, ncoeff, pdl, pd, pdu)
  
   !=======SAVE RESULTS=======

  !dgtlmpep ----------
  WRITE(1,'(20G15.4)', advance='no') timeStats(1)
  WRITE(1, '(A)', advance='no') ', '
  !------------------

  !qep3d -----------
  WRITE(1, '(20G15.4)', advance='no') timeStats(2)
  !------------------
  WRITE(1, *)
END DO
 CLOSE(UNIT=1)

END PROGRAM main




