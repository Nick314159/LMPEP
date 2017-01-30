PROGRAM sample_tri_test_driver2
USE environment
USE dgtlmpep_subroutines
USE qep3dea
USE qep3deacx
USE qep3dlag
IMPLICIT NONE

!===Variables===
!LMPEP
INTEGER :: d
INTEGER, DIMENSION(4) :: iseed
REAL(dp), DIMENSION(:), ALLOCATABLE :: cond, er, ei, ncoeff, berr, ferr
REAL(dp), DIMENSION(:,:), ALLOCATABLE :: pd, pdl, pdu, xr, xi, yr, yi
!QEP3D
INTEGER :: mode, n, i, neg, detsgn, mxit, iter, itermx, imax, jmax, jmin
INTEGER, PARAMETER :: ATTEMPTS=200
REAL(dp) :: alpha
REAL(dp), ALLOCATABLE, DIMENSION(:) :: a, b, c, au, bu, cu, al, bl, cl, z
COMPLEX(dp), ALLOCATABLE, DIMENSION(:) :: zcx, co, si, ad, adl, adu, x, y
!testing
INTEGER :: clock, clock_rate, clock_start, clock_stop, k, ppos
CHARACTER(LEN=64), DIMENSION(7) :: tests
REAL(dp), DIMENSION(2) :: timeStats
!intrinsic procedures
INTRINSIC :: COUNT, DBLE, MAXVAL, MOD, SYSTEM_CLOCK
!external procedures
REAL(dp) :: dlangt
EXTERNAL :: dlangt

!create iseed, used in dgtlmpep
CALL SYSTEM_CLOCK(COUNT=clock)
CALL srand(clock)
DO i=1,4
  iseed(i)=MOD(irand(),4095)
ENDDO
IF(MOD(iseed(4),2)==0) THEN
  iseed(4)=iseed(4)+1
ENDIF

tests(1) = 'data_Ex102_100_EAC.dat'
tests(2) = 'data_Ex102_100_EAR.dat'
tests(3) = 'data_Ex102_100_LAG.dat'
tests(4) = 'data_Ex103A_100_EAC.dat'
tests(5) = 'data_Ex103A_200_EAC.dat'
tests(6) = 'data_Ex103B_100_EAC.dat'
tests(7) = 'data_Ex103B_200_EAC.dat'

!degree parameter
d=2

!===DGTLMPEP/QEP3D Accuracy/Timing Tests===
OPEN(UNIT=1,FILE=resultsDir//"outputSampleTri.csv")
WRITE(1, '(A)',  advance='no') 'Problem,        '
WRITE(1, '(A)',  advance='no') 'DGTLMPEP TIME,   '
WRITE(1, '(A)',  advance='no') 'DGTLMPEP MAX FERR,   '
WRITE(1, '(A)',  advance='no') 'QEP3D TIME,      '
WRITE(1, '(A)',  advance='no') 'QEP3D MAX FERR,   '
WRITE(1, *)
DO k=1,7
!!! Open file
  OPEN(UNIT=2,FILE=sampleProblemsDir//tests(k))
  ppos = scan(trim(tests(k)),".", BACK= .true.) - 1
  WRITE(1, '(A)', advance='no') tests(k)(1:ppos)
  WRITE(*,*) 'Testing '//tests(k)(1:ppos)
  WRITE(1, '(A)', advance='no') ','
  READ(2,'(I2)') mode
  READ(2, '(I10)') n
!!! Allocate
  ALLOCATE(a(n),b(n),c(n),au(n-1),al(n-1),bu(n-1),bl(n-1),cu(n-1),cl(n-1))
  ALLOCATE(pdl(n-1,d+1),pd(n,d+1),pdu(n-1,d+1))
  ALLOCATE(ncoeff(d+1))
  ALLOCATE(xr(n,n*d), xi(n,n*d), yr(n,n*d), yi(n,n*d))
  ALLOCATE(berr(n*d), cond(n*d), ferr(n*d), er(n*d), ei(n*d))
  ALLOCATE(z(2*n),zcx(2*n))
  ALLOCATE(ad(n), adu(n-1), adl(n-1), co(n-1), si(n-1), x(n), y(n))
!!! Read file
  IF(mode<5) THEN
    READ(2,*) (a(i), i=1,n)
    READ(2,*) (au(i), i=1,n-1)
    READ(2,*) (b(i), i=1,n)
    READ(2,*) (bu(i), i=1,n-1)
    READ(2,*) (c(i), i=1,n)
    READ(2,*) (cu(i), i=1,n-1)
    pd(:,3)=a
    pdu(:,3)=au
    pdl(:,3)=au
    pd(:,2)=b
    pdu(:,2)=bu
    pdl(:,2)=bu
    pd(:,1)=c
    pdu(:,1)=cu
    pdl(:,1)=cu
  ELSE
    READ(2,*) (a(i), i=1,n)
    READ(2,*) (au(i), i=1,n-1)
    READ(2,*) (al(i), i=1,n-1)
    READ(2,*) (b(i), i=1,n)
    READ(2,*) (bu(i), i=1,n-1)
    READ(2,*) (bl(i), i=1,n-1)
    READ(2,*) (c(i), i=1,n)
    READ(2,*) (cu(i), i=1,n-1)
    READ(2,*) (cl(i), i=1,n-1)
    pd(:,3)=a
    pdu(:,3)=au
    pdl(:,3)=al
    pd(:,2)=b
    pdu(:,2)=bu
    pdl(:,2)=bl
    pd(:,1)=c
    pdu(:,1)=cu
    pdl(:,1)=cl
  ENDIF
  CLOSE(UNIT=2)

!!! Compute eigenvalues using DGTLMPEP
  !coefficient norm
  DO i=1,d+1
    ncoeff(i)=dlangt('F',n,pdl(1,i),pd(1,i),pdu(1,i))
  ENDDO
  !solve problem
  CALL SYSTEM_CLOCK(count_rate=clock_rate)
  CALL SYSTEM_CLOCK(COUNT=clock_start)
  CALL dgtlm(pdl, pd, pdu, xr, xi, yr, yi, er, ei, berr, ncoeff, iseed, d, n)
  CALL SYSTEM_CLOCK(COUNT=clock_stop)
  WRITE(1,'(20G15.4)', advance='no')  DBLE(clock_stop-clock_start)/DBLE(clock_rate)
  WRITE(1, '(A)', advance='no') ', '
  !error estimates for dgtlmpep
  CALL dposterrcond(pdl,pd,pdu,xr,xi,yr,yi,er,ei,ncoeff,berr,cond,ferr,d,n)
  WRITE(1,'(20G15.4)', advance='no')  MAXVAL(ferr)
  WRITE(1, '(A)', advance='no') ', '

!!! Compute eigenvalues using QEP3D
  z = zero
  zcx = zero	
  mxit = 400 ! maximal number of iteration
  iter = 0
  itermx = 500
  CALL SYSTEM_CLOCK(count_rate=clock_rate)
  CALL SYSTEM_CLOCK(COUNT=clock_start)	
  ! Real Ehrlich-Aberth
  IF (mode>=1 .AND. mode<=3) THEN
    CALL reigen(a,au,b,bu,c,cu,n,z,mxit,iter,itermx,imax,mode)
  ENDIF
  ! Laguerre
  IF (mode==4) THEN
    CALL reigenl(a,au,b,bu,c,cu,n,z,mxit,iter,itermx,imax)
  ENDIF
  ! Complex Ehrlich Aberth
  IF (mode>=5 .AND. mode<=7) THEN
    CALL reigencx(a,au,al,b,bu,bl,c,cu,cl,n,zcx,mxit,iter,itermx,imax,mode-4)
  ENDIF
  CALL SYSTEM_CLOCK(COUNT=clock_stop)
  WRITE(1,'(20G15.4)', advance='no') DBLE(clock_stop-clock_start)/DBLE(clock_rate)
  WRITE(1, '(A)', advance='no') ', '
  !eigenvectors
  IF(mode<5) zcx=DCMPLX(z)
  DO i=1,n*d
    er(i)=DBLE(zcx(i)); ei(i)=DIMAG(zcx(i))
    IF(ZABS(zcx(i))>one) THEN
      zcx(i)=1/zcx(i)
      CALL zrevgteval(pdl, pd, pdu, zcx(i), adl, ad, adu, d, n, 0)
      CALL drevseval(ncoeff, ZABS(zcx(i)), alpha, d, 0)
      adl=adl/alpha; ad=ad/alpha; adu=adu/alpha
      CALL zgtqr(adl, ad, adu, co, si, n)
      jmax=ZGTJMAX(ad,n)
      jmin=ZGTJMIN(ad,n)
      IF(ZABS(ad(jmin))<eps) THEN
        CALL zgtker1(adl, ad, adu, x, y, co, si, jmin, jmax, n)
        xr(1:n,i)=DBLE(x); xi(1:n,i)=DIMAG(x)
        yr(1:n,i)=DBLE(y); yi(1:n,i)=DIMAG(y)
      ELSE
        CALL zgtker2(adl, ad, adu, x, y, co, si, jmin, jmax, n)
        xr(1:n,i)=DBLE(x); xi(1:n,i)=DIMAG(x)
        yr(1:n,i)=DBLE(y); yi(1:n,i)=DIMAG(y)
      ENDIF
    ELSE
      CALL zgteval(pdl, pd, pdu, zcx(i), adl, ad, adu, d, n, 0)
      CALL dseval(ncoeff, ZABS(zcx(i)), alpha, d, 0)
      adl=adl/alpha; ad=ad/alpha; adu=adu/alpha
      CALL zgtqr(adl, ad, adu, co, si, n)
      jmax=ZGTJMAX(ad,n)
      jmin=ZGTJMIN(ad,n)
      IF(ZABS(ad(jmin))<eps) THEN
        CALL zgtker1(adl, ad, adu, x, y, co, si, jmin, jmax, n)
        xr(1:n,i)=DBLE(x); xi(1:n,i)=DIMAG(x)
        yr(1:n,i)=DBLE(y); yi(1:n,i)=DIMAG(y)
      ELSE
        CALL zgtker2(adl, ad, adu, x, y, co, si, jmin, jmax, n)
        xr(1:n,i)=DBLE(x); xi(1:n,i)=DIMAG(x)
        yr(1:n,i)=DBLE(y); yi(1:n,i)=DIMAG(y)
      ENDIF
    ENDIF
  ENDDO
  !error estimates for QEP3D
  CALL dposterrcond(pdl,pd,pdu,xr,xi,yr,yi,er,ei,ncoeff,berr,cond,ferr,d,n)
  WRITE(1,'(20G15.4)', advance='no')  MAXVAL(ferr)
  WRITE(1, *)

!!! Deallocate
  DEALLOCATE(a,b,c,au,al,bu,bl,cu,cl)
  DEALLOCATE(pdl,pd,pdu)
  DEALLOCATE(ncoeff)
  DEALLOCATE(xr,xi,yr,yi)
  DEALLOCATE(berr,cond,ferr,er,ei)
  DEALLOCATE(z,zcx)
  DEALLOCATE(ad, adu, adl, co, si, x, y)
ENDDO
  CLOSE(UNIT=1)


END PROGRAM sample_tri_test_driver2
