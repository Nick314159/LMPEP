PROGRAM tri_test_driver
USE environment
USE dgtlmpep_subroutines
USE qep3deacx
IMPLICIT NONE

!=======VARIABLES=======
!Common
INTEGER :: clock, clock_rate, clock_start, clock_stop, i, j, m, n, k, info
INTEGER :: d, startSize, maxSize, jumpFactor
INTEGER, DIMENSION(4) :: iseed
REAL(dp), DIMENSION(:,:), ALLOCATABLE :: timeStats
CHARACTER(LEN=100) :: arg
REAL(dp), DIMENSION(:), ALLOCATABLE :: berr1, berr2, cond, ferr, er, ei, ncoeff, work, x
REAL(dp), DIMENSION(:,:), ALLOCATABLE :: p, pdl, pd, pdu, xr, xi, yr, yi
!intrinsic procedures
INTRINSIC :: COUNT, DBLE, MAX, MAXVAL, MOD, NEW_LINE, SYSTEM_CLOCK
!external procedures
REAL(dp) :: dlangt
EXTERNAL :: dlangt
!------------------

!QEP3D
INTEGER	:: mode, neg, detsgn, jmax, jmin
INTEGER, PARAMETER :: ATTEMPTS=200
REAL(dp) :: alpha
REAL(dp), ALLOCATABLE, DIMENSION(:)	:: a, b, c, au, bu, cu, al, bl, cl
COMPLEX(dp), ALLOCATABLE, DIMENSION(:) :: z, ad, adl, adu, co, si
INTEGER	:: mxit, iter, itermx, imax
!get information
CALL GETARG(1,arg)
READ(arg, *) startSize
CALL GETARG(2,arg)
READ(arg, *) maxSize
CALL GETARG(3, arg)
READ (arg,*) m
!CALL GETARG(3,BUFFER)
!READ(BUFFER, *) jumpFactor
jumpFactor=2

ALLOCATE(timeStats(m,2))
	
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
WRITE(1, '(A)',  advance='no') 'Size,        '
WRITE(1, '(A)',  advance='no') 'DGTLMPEP TIME,   '
WRITE(1, '(A)',  advance='no') 'QEP3D TIME'
WRITE(1, *)
d = 2
n = startSize
DO WHILE (n < maxSize)
  WRITE(1, '(i6)', advance='no') n
  WRITE(1, '(A)', advance='no') ', '
  CALL SYSTEM_CLOCK(count_rate=clock_rate)
  DO j=1,m
    !create tridiagnol case
    ALLOCATE(work(n+n*(d+1)), x(n), p(n,n*(d+1)))
    DO i=1,d+1
      CALL dlarnv(2, iseed, n, x)
      CALL dlagge(n, n, 1, 1, x, p(1,n*(i-1)+1), n, iseed, work, info)
    ENDDO
    DEALLOCATE(work, x)
    
    !store tridiagonal structure
    ALLOCATE(pdl(n-1,d+1), pd(n,d+1), pdu(n-1,d+1))
    DO k=1,d+1
      !pdl
      DO i=1,n-1
        pdl(i,k)=p(i+1,(k-1)*n+i)
      ENDDO
      !pd
      DO i=1,n
        pd(i,k)=p(i,(k-1)*n+i)
      ENDDO
      !pdu
      DO i=2,n
       pdu(i-1,k)=p(i-1,(k-1)*n+i)
      ENDDO
    ENDDO
    DEALLOCATE(p)

    !solve dgtlmpep
    ALLOCATE(xr(n,n*d), xi(n,n*d), yr(n,n*d), yi(n,n*d))
    ALLOCATE(berr1(n*d), berr2(n*d), er(n*d), ei(n*d), ncoeff(d+1))
    DO i=1,d+1
      ncoeff(i)=dlangt('F',n,pdl(1,i),pd(1,i),pdu(1,i))
    ENDDO
    CALL SYSTEM_CLOCK(count_rate=clock_rate)
    CALL SYSTEM_CLOCK(COUNT=clock_start)
    CALL dgtlm(pdl, pd, pdu, xr, xi, yr, yi, er, ei, berr1, ncoeff, iseed, d, n)
    CALL SYSTEM_CLOCK(COUNT=clock_stop)
    timeStats(j,1) = DBLE(clock_stop-clock_start)/DBLE(clock_rate)
    !deallocate
    DEALLOCATE(xr, xi, yr, yi, er, ei)
    !print berr for dgtlmpep
    PRINT*, 'max berr dgtlm =', MAXVAL(berr1)
    
    ALLOCATE(a(n), b(n), c(n), au(n), bu(n), cu(n), al(n), bl(n), cl(n))
    !a
    a = pd(:, 3)
    !au
    au = pdu(:, 3)
    !al
    al = pdl(:, 3)
    !b
    b = pd(:, 2)
    !bu
    bu = pdu(:, 2)
    !bl
    bl = pdl(:, 2)
    !c
    c = pd(:, 1)
    !cu
    cu = pdu(:, 1)
    !cl
    cl = pdl(:, 1)

    !maximal number of iteration
    mxit = 400
    iter = 0
    itermx = 500
    !solve QEP3D
    ALLOCATE(z(2*n))
    CALL SYSTEM_CLOCK(count_rate=clock_rate)
    CALL SYSTEM_CLOCK(COUNT=clock_start)
    CALL reigencx(a,au,al,b,bu,bl,c,cu,cl,n,z,mxit,iter,itermx,imax,2)
    CALL SYSTEM_CLOCK(COUNT=clock_stop)
    
    timeStats(j,2) = DBLE(clock_stop-clock_start)/DBLE(clock_rate)
    DEALLOCATE(a, b, c, au, bu, cu, al, bl, cl)

    !backward error estimates for QEP3D
    ALLOCATE(ad(n), adu(n-1), adl(n-1), co(n-1), si(n-1))
    DO i=1,d*n
      IF(ZABS(z(i))>one) THEN
        CALL zrevgteval(pdl, pd, pdu, z(i), adl, ad, adu, d, n, 0)
        CALL drevseval(ncoeff, ZABS(z(i)), alpha, d, 0)
        adl=adl/alpha; ad=ad/alpha; adu=adu/alpha
        CALL zgtqr(adl, ad, adu, co, si, n)
        jmax=ZGTJMAX(ad,n)
        jmin=ZGTJMIN(ad,n)
        IF(ZABS(ad(jmin))<eps) THEN
          berr2(i)=ZABS(ad(jmin))
        ELSE
          CALL zberrapprox(adl, ad, adu, co, si, iseed, berr2(i), n)
        ENDIF
      ELSE
        CALL zgteval(pdl, pd, pdu, z(i), adl, ad, adu, d, n, 0)
        CALL dseval(ncoeff, ZABS(z(i)), alpha, d, 0)
        adl=adl/alpha; ad=ad/alpha; adu=adu/alpha
        CALL zgtqr(adl, ad, adu, co, si, n)
        jmax=ZGTJMAX(ad,n)
        jmin=ZGTJMIN(ad,n)
        IF(ZABS(ad(jmin))<eps) THEN
          berr2(i)=ZABS(ad(jmin))
        ELSE
          CALL zberrapprox(adl, ad, adu, co, si, iseed, berr2(i), n)
        ENDIF
      ENDIF
    ENDDO
    DEALLOCATE(adl, ad, adu, co, si, z)
    PRINT*, 'max berr qep3d =', MAXVAL(berr2)
    !deallocate
    DEALLOCATE(berr1, berr2, ncoeff, pdl, pd, pdu)
  ENDDO

  !=======SAVE RESULTS=======

  !dgtlmpep ----------
  WRITE(1,'(20G15.4)', advance='no') SUM(timeStats(:,1))/m
  WRITE(1, '(A)', advance='no') ', '
  !------------------

  !qep3d -----------
  WRITE(1, '(20G15.4)', advance='no') SUM(timeStats(:,2))/m
  !------------------
  WRITE(1, *)
  n = jumpFactor * n
ENDDO
CLOSE(UNIT=1)
END PROGRAM tri_test_driver
