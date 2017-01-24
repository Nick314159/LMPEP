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
REAL(dp), DIMENSION(:), ALLOCATABLE :: berr, cond, ferr, er, ei, ncoeff, work, x
REAL(dp), DIMENSION(:,:), ALLOCATABLE :: p, pdl, pd, pdu, xr, xi, yr, yi
!intrinsic procedures
INTRINSIC :: COUNT, DBLE, MAX, MAXVAL, MOD, NEW_LINE, SYSTEM_CLOCK
!external procedures
REAL(dp) :: dlangt
EXTERNAL :: dlangt
!------------------

!QEP3D
INTEGER	:: mode, neg, detsgn
INTEGER, PARAMETER :: ATTEMPTS = 200
REAL(dp), ALLOCATABLE, DIMENSION(:)	:: a, b, c, au, bu, cu, al, bl, cl
COMPLEX(dp), ALLOCATABLE, DIMENSION(:) :: z
COMPLEX(dp), ALLOCATABLE, DIMENSION(:)	:: zcx
REAL(4)	 :: T1
INTEGER	:: mxit, iter, iterlast, itermax
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
        pdl(i,k)=p(i+1,(k-1)*n+i)
      ENDDO
      !pd
      DO i=1,n
        pd(i,k)=p(i,(k-1)*n+i)
      ENDDO
      !pdu
      !DO i=2,n
       ! pdu(i-1,k)=p(i-1,(k-1)*n+i)
       pdu=pdl
      !ENDDO
    ENDDO
    DEALLOCATE(p)
    
    ALLOCATE(a(n), b(n), c(n), au(n), bu(n), cu(n), al(n), bl(n), cl(n))
    !a
    a = pd(:, 1)
    !au
    au = pdu(:, 1)
    !al
    al = pdl(:, 1)
    !b
    b = pd(:, 2)
    !bu
    bu = pdu(:, 2)
    !bl
    bl = pdl(:, 2)
    !c
    c = pd(:, 3)
    !cu
    cu = pdu(:, 3)
    !al
    cl = pdl(:, 3)

    !maximal number of iteration
    mxit = 400	
    iter = 0
    iterlast = 500
    !solve QEP3D
    ALLOCATE(z(2*n), zcx(2*n))
    CALL SYSTEM_CLOCK(count_rate=clock_rate)
    CALL SYSTEM_CLOCK(COUNT=clock_start)
    CALL reigencx(a,au,al,b,bu,bl,c,cu,cl,n,z,mxit,iter,iterlast,itermax,mode-4)
    CALL SYSTEM_CLOCK(COUNT=clock_stop)
    DEALLOCATE(z,zcx)
    
    timeStats(j,2) = DBLE(clock_stop-clock_start)/DBLE(clock_rate)
    DEALLOCATE(a, b, c, au, bu, cu, al, bl, cl)
    
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

    timeStats(j,1) = DBLE(clock_stop-clock_start)/DBLE(clock_rate)
 
   
    !compute backward error and condition number
    CALL dposterrcond(pdl, pd, pdu, xr, xi, yr, yi, er, ei, ncoeff, berr, cond, ferr, d, n)
    
    !deallocate
    DEALLOCATE(xr, xi, yr, yi)
    DEALLOCATE(berr, cond, ferr, er, ei, ncoeff, pdl, pd, pdu)
  
  END DO

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
END DO
CLOSE(UNIT=1)
END PROGRAM tri_test_driver
