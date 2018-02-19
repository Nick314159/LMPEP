!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!   Thomas R. Cameron* and Nikolas I. Steckley#
!                                                               
!   *Dept. Mathematics and Computer Science, Davidson College
!   #Dept. DiscoverOrg LLC.
!       
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!   Last modified 1 January 2018
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!   Test Comparison of DSLM against POLZEROS and AMVW
!   
!   Test 1: root approximations are computed for random polynomials
!   of degree startDegree to endDegree over several iterations.
!   The average computation time and backward error is recorded.
!   
!   Test 2: root approximations are computed for several famous
!   polynomials and the 2-norm of the relative forward error is recorded.
!   
!   The results for Test1 are written to
!   file: "results/testComparison1.csv".
!   
!   The results for Test2 are written to
!   file: "results/testComparison2.csv".
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!   Inpute Parameter
!   
!   1) Start Degree, default 100
!   
!   2) End Degree, default 6400
!   
!   3) Iterations, default 10
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
PROGRAM testComparison
USE poly_zeroes
IMPLICIT NONE
INTEGER                                         :: clock, clock_rate, clock_start, clock_stop
INTEGER                                         :: i, j, nitmax, iter
DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE     :: radius
DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE   :: time, error
DOUBLE PRECISION                                :: eps, big, small
DOUBLE COMPLEX, DIMENSION(:), ALLOCATABLE       :: root, poly, exacteigs
LOGICAL, DIMENSION(:), ALLOCATABLE              :: err
INTEGER                                         :: it, itmax
INTEGER                                         :: deg, startDegree, endDegree, flag
DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE     :: berr, er, ei, p, p2, alpha
CHARACTER(LEN=32)                               :: arg
DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE     :: REIGS, IEIGS
INTEGER, DIMENSION(:), ALLOCATABLE              :: ITS
DOUBLE PRECISION                                :: t2
DOUBLE COMPLEX                                  :: a, t
! Parameter
DOUBLE PRECISION, PARAMETER         						:: pi = 3.141592653589793D0
! External subroutines
DOUBLE PRECISION                                :: dznrm2
EXTERNAL                                        :: dznrm2

CALL init_random_seed()

eps    = epsilon(1.0D0)
small  = tiny(1.0D0)
big    = huge(1.0D0)
nitmax = 60
FLAG 	 = 1

!****************************************************
!                   Test 1                          *
!****************************************************
! Set Start and End Degree, and Iterations
    CALL GET_COMMAND_ARGUMENT(1,arg,status=flag)
    IF(flag==0) THEN
        READ(arg, '(I10)') startDegree
    ELSE
        startDegree=100
    ENDIF
    CALL GET_COMMAND_ARGUMENT(2,arg,status=flag)
    IF(flag==0) THEN
        READ(arg, '(I10)') endDegree
    ELSE
        endDegree=6400
    ENDIF
    CALL GET_COMMAND_ARGUMENT(3,arg,status=flag)
    IF(flag==0) THEN
        READ(arg, '(I10)') itmax
    ELSE
        itmax=10
    ENDIF

deg=startDegree
OPEN(UNIT=1,FILE="results/testComparison1.csv")
WRITE(1,'(A)') 'Degree, DSLM Time, DSLM Berr, Pzeros Time, Pzeros Berr, AMVW Time, AMVW Berr'
ALLOCATE(time(itmax, 3), error(itmax, 3))
DO WHILE(deg<=endDegree)
  WRITE(1,'(I10)', advance='no') deg
  WRITE(1,'(A)', advance='no') ','
  
  DO it = 1, itmax
		!Poly, coeff moduli, and berr
		ALLOCATE(p(deg+1), berr(deg), alpha(deg+1))
    CALL daruv(deg+1,p)
    DO j = 1, deg+1
       alpha(j)=dabs(p(j))*(4*j-3)
    END DO
    !LMPEP
    ALLOCATE(er(deg), ei(deg))
    CALL system_clock(count_rate=clock_rate)
    CALL system_clock(count=clock_start)
    CALL dslm(p, deg, er, ei, berr)
    CALL system_clock(count=clock_stop)
    time(it, 1)=(dble(clock_stop-clock_start)/dble(clock_rate))
    error(it, 1) = maxval(berr)
		!Deallocate LMPEP
		DEALLOCATE(er, ei)

    !Pzeros
    ALLOCATE(poly(0:deg), radius(1:deg), root(1:deg), err(deg+1)) 
    DO i= 0, deg
      poly(i)=cmplx(p(i+1), 0.0D0, dp)
    ENDDO
    CALL system_clock(count_rate=clock_rate)
    CALL system_clock(count=clock_start)
    CALL polzeros(deg, poly, eps, big, small, nitmax, root, radius, err, iter)
    CALL system_clock(count=clock_stop)
    time(it, 2)=(dble(clock_stop-clock_start)/dble(clock_rate))
    DO j = 1, deg
      t =  root(j)
      t2 = zabs(t)
      IF (t2>1) THEN
        t = t**(-1)
        t2 = t2**(-1)
        CALL dzrevseval(p, t, deg, 0, a)
        CALL drevseval(alpha, t2, deg, 0, berr(j))
      ELSE
        CALL dzseval(p, t, deg, 0, a)
        CALL dseval(alpha, t2, deg, 0, berr(j)) 
      END IF
      berr(j) = zabs(a)/berr(j)
    END DO
    error(it, 2) = maxval(berr)
		!Deallocate Pzeros
    DEALLOCATE(poly, radius ,root, err)

    !AMVW
    ALLOCATE(REIGS(deg),IEIGS(deg),ITS(deg), p2(deg))
    DO i=1,deg
        p2(deg-i+1)=p(i)/p(deg+1)
    ENDDO
    CALL system_clock(count_rate=clock_rate)
    CALL system_clock(count=clock_start)
    CALL damvw(deg, p2, REIGS, IEIGS, ITS, FLAG)
    CALL system_clock(count=clock_stop)  
    time(it, 3) = dble(clock_stop-clock_start)/dble(clock_rate)
    DO j = 1, deg
      t =  cmplx(REIGS(j), IEIGS(j),dp)
      t2 = zabs(t)
      IF (t2>1) THEN
        t = t**(-1)
        t2 = t2**(-1)
        CALL dzrevseval(p, t, deg, 0, a)
        CALL drevseval(alpha, t2, deg, 0, berr(j))
      ELSE
        CALL dzseval(p, t, deg, 0, a)
        CALL dseval(alpha, t2, deg, 0, berr(j)) 
      END IF
      berr(j) = zabs(a)/berr(j)
    END DO
    error(it, 3) = maxval(berr)
		!Deallocate AMVW
    DEALLOCATE(REIGS, IEIGS, ITS, p2)
		!Deallocate poly
		DEALLOCATE(p, berr, alpha)
  ENDDO
  
  WRITE(1,'(ES15.2)', advance='no') sum(time(:,1))/itmax
  WRITE(1,'(A)', advance='no') ','
  WRITE(1,'(ES15.2)', advance='no') sum(error(:,1))/itmax
  WRITE(1,'(A)', advance='no') ','
  WRITE(1,'(ES15.2)', advance='no') sum(time(:,2))/itmax
  WRITE(1,'(A)', advance='no') ',' 
  WRITE(1,'(ES15.2)', advance='no') sum(error(:,2))/itmax
  WRITE(1,'(A)', advance='no') ','
  WRITE(1,'(ES15.2)', advance='no') sum(time(:,3))/itmax
  WRITE(1,'(A)', advance='no') ','
  WRITE(1,'(ES15.2)') sum(error(:,3))/itmax
  deg=2*deg
ENDDO
DEALLOCATE(time, error)
    CLOSE(1)

!****************************************************
!                   Test 2                          *
!****************************************************
! Open results file
OPEN(UNIT=1,FILE="results/testComparison2.csv")
WRITE(1,'(A)') ' , DSLM, POLZEROS, AMVW'

! Test 1: Wilkinson Deg 10 Polynomial
deg=10
ALLOCATE(exacteigs(deg), p(deg+1))
p(1)=3628800D0
p(2)=-10628640D0
p(3)=12753576D0
p(4)=-8409500D0
p(5)=3416930D0
p(6)=-902055D0
p(7)=157773D0
p(8)=-18150D0
p(9)=1320D0
p(10)=-55D0
p(11)=1D0
exacteigs(1)=cmplx(1D0,0D0,dp)
exacteigs(2)=cmplx(2D0,0D0,dp)
exacteigs(3)=cmplx(3D0,0D0,dp)
exacteigs(4)=cmplx(4D0,0D0,dp)
exacteigs(5)=cmplx(5D0,0D0,dp)
exacteigs(6)=cmplx(6D0,0D0,dp)
exacteigs(7)=cmplx(7D0,0D0,dp)
exacteigs(8)=cmplx(8D0,0D0,dp)
exacteigs(9)=cmplx(9D0,0D0,dp)
exacteigs(10)=cmplx(10D0,0D0,dp)
WRITE(1, '(A)', advance='no')  'Test 1, '
! DSLM
ALLOCATE(er(deg),ei(deg),berr(deg))
CALL dslm(p, deg, er, ei, berr)
CALL dsort(er, ei, deg)
WRITE(1, '(ES15.2)', advance='no') dznrm2(deg,cmplx(er,ei,dp)-exacteigs,1)/dznrm2(deg,exacteigs,1)
WRITE(1,'(A)',advance='no') ','
DEALLOCATE(er,ei,berr)
! POLZEROS
ALLOCATE(poly(0:deg), radius(1:deg), root(1:deg), err(deg+1)) 
DO i= 0, deg
  poly(i)=cmplx(p(i+1), 0.0D0, dp)
ENDDO
CALL polzeros(deg, poly, eps, big, small, nitmax, root, radius, err, iter)
CALL zsort(root, deg)
WRITE(1, '(ES15.2)', advance='no') dznrm2(deg,root-exacteigs,1)/dznrm2(deg,exacteigs,1)
WRITE(1,'(A)', advance='no') ','
DEALLOCATE(poly,radius,root,err)
! AMVW
ALLOCATE(REIGS(deg),IEIGS(deg),ITS(deg), p2(deg))
DO i=1,deg
    p2(deg-i+1)=p(i)
ENDDO
CALL damvw(deg, p2, REIGS, IEIGS, ITS, FLAG)
CALL dsort(REIGS,IEIGS,deg)
WRITE(1, '(ES15.2)') dznrm2(deg,cmplx(REIGS,IEIGS,dp)-exacteigs,1)/dznrm2(deg,exacteigs,1)
DEALLOCATE(REIGS,IEIGS,ITS,p2)
! End Test 1
DEALLOCATE(exacteigs, p)

! Test 2: Wilkinson Deg 15 Polynomial
deg=15
ALLOCATE(exacteigs(deg),p(deg+1))
p(1)=-1307674368000D0
p(2)=4339163001600D0
p(3)=-6165817614720D0
p(4)=5056995703824D0
p(5)=-2706813345600D0
p(6)=1009672107080D0
p(7)=-272803210680D0
p(8)=54631129553D0
p(9)=-8207628000D0
p(10)=928095740D0
p(11)=-78558480D0
p(12)=4899622D0
p(13)=-218400D0
p(14)=6580D0
p(15)=-120D0
p(16)=1D0
exacteigs(1)=dcmplx(1D0,0D0)
exacteigs(2)=dcmplx(2D0,0D0)
exacteigs(3)=dcmplx(3D0,0D0)
exacteigs(4)=dcmplx(4D0,0D0)
exacteigs(5)=dcmplx(5D0,0D0)
exacteigs(6)=dcmplx(6D0,0D0)
exacteigs(7)=dcmplx(7D0,0D0)
exacteigs(8)=dcmplx(8D0,0D0)
exacteigs(9)=dcmplx(9D0,0D0)
exacteigs(10)=dcmplx(10D0,0D0)
exacteigs(11)=dcmplx(11D0,0D0)
exacteigs(12)=dcmplx(12D0,0D0)
exacteigs(13)=dcmplx(13D0,0D0)
exacteigs(14)=dcmplx(14D0,0D0)
exacteigs(15)=dcmplx(15D0,0D0)
WRITE(1, '(A)', advance='no')  'Test 2, '
! DSLM
ALLOCATE(er(deg),ei(deg),berr(deg))
CALL dslm(p, deg, er, ei, berr)
CALL dsort(er, ei, deg)
WRITE(1, '(ES15.2)', advance='no') dznrm2(deg,dcmplx(er,ei)-exacteigs,1)/dznrm2(deg,exacteigs,1)
WRITE(1,'(A)',advance='no') ','
DEALLOCATE(er,ei,berr)
! POLZEROS
ALLOCATE(poly(0:deg), radius(1:deg), root(1:deg), err(deg+1)) 
DO i= 0, deg
  poly(i)=dcmplx(p(i+1), 0D0)
ENDDO
CALL polzeros(deg, poly, eps, big, small, nitmax, root, radius, err, iter)
CALL zsort(root, deg)
WRITE(1, '(ES15.2)', advance='no') dznrm2(deg,root-exacteigs,1)/dznrm2(deg,exacteigs,1)
WRITE(1,'(A)', advance='no') ','
DEALLOCATE(poly,radius,root,err)
! AMVW
ALLOCATE(REIGS(deg),IEIGS(deg),ITS(deg), p2(deg))
DO i=1,deg
    p2(deg-i+1)=p(i)
ENDDO
CALL damvw(deg, p2, REIGS, IEIGS, ITS, FLAG)
CALL dsort(REIGS,IEIGS,deg)
WRITE(1, '(ES15.2)') dznrm2(deg,dcmplx(REIGS,IEIGS)-exacteigs,1)/dznrm2(deg,exacteigs,1)
DEALLOCATE(REIGS,IEIGS,ITS,p2)
! End Test 2
DEALLOCATE(exacteigs, p)

! Test 3: Wilkinson Deg 20 Polynomial
deg=20
ALLOCATE(exacteigs(deg),p(deg+1))
p(1)=2432902008176640000D0
p(2)=-8752948036761600000D0
p(3)=13803759753640704000D0
p(4)=-12870931245150988800D0
p(5)=8037811822645051776D0
p(6)=-3599979517947607200D0
p(7)=1206647803780373360D0
p(8)=-311333643161390640D0
p(9)=63030812099294896D0
p(10)=-10142299865511450D0
p(11)=1307535010540395D0
p(12)=-135585182899530D0
p(13)=11310276995381D0
p(14)=-756111184500D0
p(15)=40171771630D0
p(16)=-1672280820D0
p(17)=53327946D0
p(18)=-1256850D0
p(19)=20615D0
p(20)=-210D0
p(21)=1D0
exacteigs(1)=dcmplx(1D0,0D0)
exacteigs(2)=dcmplx(2D0,0D0)
exacteigs(3)=dcmplx(3D0,0D0)
exacteigs(4)=dcmplx(4D0,0D0)
exacteigs(5)=dcmplx(5D0,0D0)
exacteigs(6)=dcmplx(6D0,0D0)
exacteigs(7)=dcmplx(7D0,0D0)
exacteigs(8)=dcmplx(8D0,0D0)
exacteigs(9)=dcmplx(9D0,0D0)
exacteigs(10)=dcmplx(10D0,0D0)
exacteigs(11)=dcmplx(11D0,0D0)
exacteigs(12)=dcmplx(12D0,0D0)
exacteigs(13)=dcmplx(13D0,0D0)
exacteigs(14)=dcmplx(14D0,0D0)
exacteigs(15)=dcmplx(15D0,0D0)
exacteigs(16)=dcmplx(16D0,0D0)
exacteigs(17)=dcmplx(17D0,0D0)
exacteigs(18)=dcmplx(18D0,0D0)
exacteigs(19)=dcmplx(19D0,0D0)
exacteigs(20)=dcmplx(20D0,0D0)
WRITE(1, '(A)', advance='no')  'Test 3, '
! DSLM
ALLOCATE(er(deg),ei(deg),berr(deg))
CALL dslm(p, deg, er, ei, berr)
CALL dsort(er, ei, deg)
WRITE(1, '(ES15.2)', advance='no') dznrm2(deg,dcmplx(er,ei)-exacteigs,1)/dznrm2(deg,exacteigs,1)
WRITE(1,'(A)',advance='no') ','
DEALLOCATE(er,ei,berr)
! POLZEROS
ALLOCATE(poly(0:deg), radius(1:deg), root(1:deg), err(deg+1)) 
DO i= 0, deg
  poly(i)=dcmplx(p(i+1), 0D0)
ENDDO
CALL polzeros(deg, poly, eps, big, small, nitmax, root, radius, err, iter)
CALL zsort(root, deg)
WRITE(1, '(ES15.2)', advance='no') dznrm2(deg,root-exacteigs,1)/dznrm2(deg,exacteigs,1)
WRITE(1,'(A)', advance='no') ','
DEALLOCATE(poly,radius,root,err)
! AMVW
ALLOCATE(REIGS(deg),IEIGS(deg),ITS(deg), p2(deg))
DO i=1,deg
    p2(deg-i+1)=p(i)
ENDDO
CALL damvw(deg, p2, REIGS, IEIGS, ITS, FLAG)
CALL dsort(REIGS,IEIGS,deg)
WRITE(1, '(ES15.2)') dznrm2(deg,dcmplx(REIGS,IEIGS)-exacteigs,1)/dznrm2(deg,exacteigs,1)
DEALLOCATE(REIGS,IEIGS,ITS,p2)
! End Test 3
DEALLOCATE(exacteigs, p)

! Test 4: scaled and shifted Wilkinson Deg 20 Polynomial
deg=20
ALLOCATE(exacteigs(deg),p(deg+1))
p(21)=1D0
p(20)=4D0
p(19)=-57D0/10D0
p(18)=-969D0/25D0
p(17)=-10659D0/10000D0
p(16)=93993D0/625D0
p(15)=1867909D0/25000D0
p(14)=-18824117D0/62500D0
p(13)=-10834566327D0/50000000D0
p(12)=4187079039D0/12500000D0
p(11)=141135792837D0/500000000D0
p(10)=-259286608191D0/1250000000D0
p(9)=-94957438655047D0/500000000000D0
p(8)=4248959006581D0/62500000000D0
p(7)=163007395518693D0/2500000000000D0
p(6)=-65260688438889D0/6250000000000D0
p(5)=-102724049585427219D0/10000000000000000D0
p(4)=11456073304317D0/20000000000000D0
p(3)=2280736816325919D0/4000000000000000D0
p(2)=-1899923154129D0/400000000000000D0
p(1)=-758069338497471D0/160000000000000000D0
exacteigs(1)=-2.1D0
exacteigs(2)=-1.9D0
exacteigs(3)=-1.7D0
exacteigs(4)=-1.5D0
exacteigs(5)=-1.3D0
exacteigs(6)=-1.1D0
exacteigs(7)=-0.9D0
exacteigs(8)=-0.7D0
exacteigs(9)=-0.5D0
exacteigs(10)=-0.3D0
exacteigs(11)=-0.1D0
exacteigs(12)=0.1D0
exacteigs(13)=0.3D0
exacteigs(14)=0.5D0
exacteigs(15)=0.7D0
exacteigs(16)=0.9D0
exacteigs(17)=1.1D0
exacteigs(18)=1.3D0
exacteigs(19)=1.5D0
exacteigs(20)=1.7D0
WRITE(1, '(A)', advance='no')  'Test 4, '
! DSLM
ALLOCATE(er(deg),ei(deg),berr(deg))
CALL dslm(p, deg, er, ei, berr)
CALL dsort(er, ei, deg)
WRITE(1, '(ES15.2)', advance='no') dznrm2(deg,dcmplx(er,ei)-exacteigs,1)/dznrm2(deg,exacteigs,1)
WRITE(1,'(A)',advance='no') ','
DEALLOCATE(er,ei,berr)
! POLZEROS
ALLOCATE(poly(0:deg), radius(1:deg), root(1:deg), err(deg+1)) 
DO i= 0, deg
  poly(i)=dcmplx(p(i+1), 0D0)
ENDDO
CALL polzeros(deg, poly, eps, big, small, nitmax, root, radius, err, iter)
CALL zsort(root, deg)
WRITE(1, '(ES15.2)', advance='no') dznrm2(deg,root-exacteigs,1)/dznrm2(deg,exacteigs,1)
WRITE(1,'(A)', advance='no') ','
DEALLOCATE(poly,radius,root,err)
! AMVW
ALLOCATE(REIGS(deg),IEIGS(deg),ITS(deg), p2(deg))
DO i=1,deg
    p2(deg-i+1)=p(i)
ENDDO
CALL damvw(deg, p2, REIGS, IEIGS, ITS, FLAG)
CALL dsort(REIGS,IEIGS,deg)
WRITE(1, '(ES15.2)') dznrm2(deg,dcmplx(REIGS,IEIGS)-exacteigs,1)/dznrm2(deg,exacteigs,1)
DEALLOCATE(REIGS,IEIGS,ITS,p2)
! End Test 4
DEALLOCATE(exacteigs, p)


! Test 5: reverse Wilkinson Deg 10 polynomial
deg=10
ALLOCATE(exacteigs(deg),p(deg+1))
p(11)=1D0
p(10)=-7381D0/2520D0
p(9)=177133D0/50400D0
p(8)=-84095D0/36288D0
p(7)=341693D0/362880D0
p(6)=-8591D0/34560D0
p(5)=7513D0/172800D0
p(4)=-121D0/24192D0
p(3)=11D0/30240D0
p(2)=-11D0/725760D0
p(1)=1D0/3628800D0
exacteigs(10)=dcmplx(1D0,0D0)
exacteigs(9)=dcmplx(1D0/2D0,0D0)
exacteigs(8)=dcmplx(1D0/3D0,0D0)
exacteigs(7)=dcmplx(1D0/4D0,0D0)
exacteigs(6)=dcmplx(1D0/5D0,0D0)
exacteigs(5)=dcmplx(1D0/6D0,0D0)
exacteigs(4)=dcmplx(1D0/7D0,0D0)
exacteigs(3)=dcmplx(1D0/8D0,0D0)
exacteigs(2)=dcmplx(1D0/9D0,0D0)
exacteigs(1)=dcmplx(1D0/10D0,0D0)
WRITE(1, '(A)', advance='no')  'Test 5, '
! DSLM
ALLOCATE(er(deg),ei(deg),berr(deg))
CALL dslm(p, deg, er, ei, berr)
CALL dsort(er, ei, deg)
WRITE(1, '(ES15.2)', advance='no') dznrm2(deg,dcmplx(er,ei)-exacteigs,1)/dznrm2(deg,exacteigs,1)
WRITE(1,'(A)',advance='no') ','
DEALLOCATE(er,ei,berr)
! POLZEROS
ALLOCATE(poly(0:deg), radius(1:deg), root(1:deg), err(deg+1)) 
DO i= 0, deg
  poly(i)=dcmplx(p(i+1), 0D0)
ENDDO
CALL polzeros(deg, poly, eps, big, small, nitmax, root, radius, err, iter)
CALL zsort(root, deg)
WRITE(1, '(ES15.2)', advance='no') dznrm2(deg,root-exacteigs,1)/dznrm2(deg,exacteigs,1)
WRITE(1,'(A)', advance='no') ','
DEALLOCATE(poly,radius,root,err)
! AMVW
ALLOCATE(REIGS(deg),IEIGS(deg),ITS(deg), p2(deg))
DO i=1,deg
    p2(deg-i+1)=p(i)
ENDDO
CALL damvw(deg, p2, REIGS, IEIGS, ITS, FLAG)
CALL dsort(REIGS,IEIGS,deg)
WRITE(1, '(ES15.2)') dznrm2(deg,dcmplx(REIGS,IEIGS)-exacteigs,1)/dznrm2(deg,exacteigs,1)
DEALLOCATE(REIGS,IEIGS,ITS,p2)
! End Test 5
DEALLOCATE(exacteigs, p)

! Test 6: reverse Wilkinson Deg 15 polynomial
deg=15
ALLOCATE(exacteigs(deg),p(deg+1))
p(16)=1D0
p(15)=-1195757D0/360360D0
p(14)=13215487D0/2802800D0
p(13)=-35118025721D0/9081072000D0
p(12)=2065639D0/997920D0
p(11)=-277382447D0/359251200D0
p(10)=2271089D0/10886400D0
p(9)=-54576553D0/1306368000D0
p(8)=4783D0/762048D0
p(7)=-324509D0/457228800D0
p(6)=109D0/1814400D0
p(5)=-26921D0/7185024000D0
p(4)=1D0/5987520D0
p(3)=-47D0/9340531200D0
p(2)=1D0/10897286400D0
p(1)=-1D0/1307674368000D0
exacteigs(15)=dcmplx(1D0,0D0)
exacteigs(14)=dcmplx(1D0/2D0,0D0)
exacteigs(13)=dcmplx(1D0/3D0,0D0)
exacteigs(12)=dcmplx(1D0/4D0,0D0)
exacteigs(11)=dcmplx(1D0/5D0,0D0)
exacteigs(10)=dcmplx(1D0/6D0,0D0)
exacteigs(9)=dcmplx(1D0/7D0,0D0)
exacteigs(8)=dcmplx(1D0/8D0,0D0)
exacteigs(7)=dcmplx(1D0/9D0,0D0)
exacteigs(6)=dcmplx(1D0/10D0,0D0)
exacteigs(5)=dcmplx(1D0/11D0,0D0)
exacteigs(4)=dcmplx(1D0/12D0,0D0)
exacteigs(3)=dcmplx(1D0/13D0,0D0)
exacteigs(2)=dcmplx(1D0/14D0,0D0)
exacteigs(1)=dcmplx(1D0/15D0,0D0)
WRITE(1, '(A)', advance='no')  'Test 6, '
! DSLM
ALLOCATE(er(deg),ei(deg),berr(deg))
CALL dslm(p, deg, er, ei, berr)
CALL dsort(er, ei, deg)
WRITE(1, '(ES15.2)', advance='no') dznrm2(deg,dcmplx(er,ei)-exacteigs,1)/dznrm2(deg,exacteigs,1)
WRITE(1,'(A)',advance='no') ','
DEALLOCATE(er,ei,berr)
! POLZEROS
ALLOCATE(poly(0:deg), radius(1:deg), root(1:deg), err(deg+1)) 
DO i= 0, deg
  poly(i)=dcmplx(p(i+1), 0D0)
ENDDO
CALL polzeros(deg, poly, eps, big, small, nitmax, root, radius, err, iter)
CALL zsort(root, deg)
WRITE(1, '(ES15.2)', advance='no') dznrm2(deg,root-exacteigs,1)/dznrm2(deg,exacteigs,1)
WRITE(1,'(A)', advance='no') ','
DEALLOCATE(poly,radius,root,err)
! AMVW
ALLOCATE(REIGS(deg),IEIGS(deg),ITS(deg), p2(deg))
DO i=1,deg
    p2(deg-i+1)=p(i)
ENDDO
CALL damvw(deg, p2, REIGS, IEIGS, ITS, FLAG)
CALL dsort(REIGS,IEIGS,deg)
WRITE(1, '(ES15.2)') dznrm2(deg,dcmplx(REIGS,IEIGS)-exacteigs,1)/dznrm2(deg,exacteigs,1)
DEALLOCATE(REIGS,IEIGS,ITS,p2)
! End Test 6
DEALLOCATE(exacteigs, p)

! Test 7: reverse Wilkinson Deg 20 polynomial
deg=20
ALLOCATE(exacteigs(deg), p(deg+1))
p(21)=1D0
p(20)=-55835135D0/15519504D0
p(19)=665690574539D0/117327450240D0
p(18)=-13334148911D0/2520460800D0
p(17)=52460655692911D0/15878903040000D0
p(16)=-31591404263D0/21349785600D0
p(15)=6670985204447D0/13450364928000D0
p(14)=-573738838201D0/4483454976000D0
p(13)=12196364570297D0/470762772480000D0
p(12)=-109542331D0/26276659200D0
p(11)=2965638101D0/5518098432000D0
p(10)=-7321967D0/131383296000D0
p(9)=384794917D0/82771476480000D0
p(8)=-31849D0/102478970880D0
p(7)=12437081D0/753220435968000D0
p(6)=-587D0/853991424000D0
p(5)=3931D0/179338199040000D0
p(4)=-1D0/1935713894400D0
p(3)=31D0/3658499260416000D0
p(2)=-1D0/11585247657984000D0
p(1)=1D0/2432902008176640000D0
exacteigs(20)=dcmplx(1D0,0D0)
exacteigs(19)=dcmplx(1D0/2D0,0D0)
exacteigs(18)=dcmplx(1D0/3D0,0D0)
exacteigs(17)=dcmplx(1D0/4D0,0D0)
exacteigs(16)=dcmplx(1D0/5D0,0D0)
exacteigs(15)=dcmplx(1D0/6D0,0D0)
exacteigs(14)=dcmplx(1D0/7D0,0D0)
exacteigs(13)=dcmplx(1D0/8D0,0D0)
exacteigs(12)=dcmplx(1D0/9D0,0D0)
exacteigs(11)=dcmplx(1D0/10D0,0D0)
exacteigs(10)=dcmplx(1D0/11D0,0D0)
exacteigs(9)=dcmplx(1D0/12D0,0D0)
exacteigs(8)=dcmplx(1D0/13D0,0D0)
exacteigs(7)=dcmplx(1D0/14D0,0D0)
exacteigs(6)=dcmplx(1D0/15D0,0D0)
exacteigs(5)=dcmplx(1D0/16D0,0D0)
exacteigs(4)=dcmplx(1D0/17D0,0D0)
exacteigs(3)=dcmplx(1D0/18D0,0D0)
exacteigs(2)=dcmplx(1D0/19D0,0D0)
exacteigs(1)=dcmplx(1D0/20D0,0D0)
WRITE(1, '(A)', advance='no')  'Test 7, '
! DSLM
ALLOCATE(er(deg),ei(deg),berr(deg))
CALL dslm(p, deg, er, ei, berr)
CALL dsort(er, ei, deg)
WRITE(1, '(ES15.2)', advance='no') dznrm2(deg,dcmplx(er,ei)-exacteigs,1)/dznrm2(deg,exacteigs,1)
WRITE(1,'(A)',advance='no') ','
DEALLOCATE(er,ei,berr)
! POLZEROS
ALLOCATE(poly(0:deg), radius(1:deg), root(1:deg), err(deg+1)) 
DO i= 0, deg
  poly(i)=dcmplx(p(i+1), 0D0)
ENDDO
CALL polzeros(deg, poly, eps, big, small, nitmax, root, radius, err, iter)
CALL zsort(root, deg)
WRITE(1, '(ES15.2)', advance='no') dznrm2(deg,root-exacteigs,1)/dznrm2(deg,exacteigs,1)
WRITE(1,'(A)', advance='no') ','
DEALLOCATE(poly,radius,root,err)
! AMVW
ALLOCATE(REIGS(deg),IEIGS(deg),ITS(deg), p2(deg))
DO i=1,deg
    p2(deg-i+1)=p(i)
ENDDO
CALL damvw(deg, p2, REIGS, IEIGS, ITS, FLAG)
CALL dsort(REIGS,IEIGS,deg)
WRITE(1, '(ES15.2)') dznrm2(deg,dcmplx(REIGS,IEIGS)-exacteigs,1)/dznrm2(deg,exacteigs,1)
DEALLOCATE(REIGS,IEIGS,ITS,p2)
! End Test 7
DEALLOCATE(exacteigs, p)

! Test 8: prescribed roots of varying scale
deg=20
ALLOCATE(exacteigs(deg), p(deg+1))
p(21)=1D0
p(20)=-1048575D0/1024D0
p(19)=183251413675D0/524288D0
p(18)=-6862582190715075D0/134217728D0
p(17)=59965700687947706355D0/17179869184D0
p(16)=-126769425631762997934675D0/1099511627776D0
p(15)=65934186820253621481357075D0/35184372088832D0
p(14)=-8505510099812717171095062675D0/562949953421312D0
p(13)=273210326382611632738979052435D0/4503599627370496D0
p(12)=-2189425218271613769209626653075D0/18014398509481984D0
p(11)=4380990637147598617372537398675D0/36028797018963968D0
p(10)=-2189425218271613769209626653075D0/36028797018963968D0
p(9)=273210326382611632738979052435D0/18014398509481984D0
p(8)=-8505510099812717171095062675D0/4503599627370496D0
p(7)=65934186820253621481357075D0/562949953421312D0
p(6)=-126769425631762997934675D0/35184372088832D0
p(5)=59965700687947706355D0/1099511627776D0
p(4)=-6862582190715075D0/17179869184D0
p(3)=183251413675D0/134217728D0
p(2)=-1048575D0/524288D0
p(1)=1D0/1024D0
exacteigs(1)=2**(-10)
exacteigs(2)=2**(-9)
exacteigs(3)=2**(-8)
exacteigs(4)=2**(-7)
exacteigs(5)=2**(-6)
exacteigs(6)=2**(-5)
exacteigs(7)=2**(-4)
exacteigs(8)=2**(-3)
exacteigs(9)=2**(-2)
exacteigs(10)=2**(-1)
exacteigs(11)=2**(0)
exacteigs(12)=2**(1)
exacteigs(13)=2**(2)
exacteigs(14)=2**(3)
exacteigs(15)=2**(4)
exacteigs(16)=2**(5)
exacteigs(17)=2**(6)
exacteigs(18)=2**(7)
exacteigs(19)=2**(8)
exacteigs(20)=2**(9)
WRITE(1, '(A)', advance='no')  'Test 8, '
! DSLM
ALLOCATE(er(deg),ei(deg),berr(deg))
CALL dslm(p, deg, er, ei, berr)
CALL dsort(er, ei, deg)
WRITE(1, '(ES15.2)', advance='no') dznrm2(deg,dcmplx(er,ei)-exacteigs,1)/dznrm2(deg,exacteigs,1)
WRITE(1,'(A)',advance='no') ','
DEALLOCATE(er,ei,berr)
! POLZEROS
ALLOCATE(poly(0:deg), radius(1:deg), root(1:deg), err(deg+1)) 
DO i= 0, deg
  poly(i)=dcmplx(p(i+1), 0D0)
ENDDO
CALL polzeros(deg, poly, eps, big, small, nitmax, root, radius, err, iter)
CALL zsort(root, deg)
WRITE(1, '(ES15.2)', advance='no') dznrm2(deg,root-exacteigs,1)/dznrm2(deg,exacteigs,1)
WRITE(1,'(A)', advance='no') ','
DEALLOCATE(poly,radius,root,err)
! AMVW
ALLOCATE(REIGS(deg),IEIGS(deg),ITS(deg), p2(deg))
DO i=1,deg
    p2(deg-i+1)=p(i)
ENDDO
CALL damvw(deg, p2, REIGS, IEIGS, ITS, FLAG)
CALL dsort(REIGS,IEIGS,deg)
WRITE(1, '(ES15.2)') dznrm2(deg,dcmplx(REIGS,IEIGS)-exacteigs,1)/dznrm2(deg,exacteigs,1)
DEALLOCATE(REIGS,IEIGS,ITS,p2)
! End Test 8
DEALLOCATE(exacteigs, p)

! Test 9: prescribed roots of varying scale - 3
deg=20
ALLOCATE(exacteigs(deg), p(deg+1))
p(21)=1D0
p(20)=-987135D0/1024D0
p(19)=153546333355D0/524288D0
p(18)=-4536701774077635D0/134217728D0
p(17)=22981827635371262835D0/17179869184D0
p(16)=-3248268654710998505043D0/1099511627776D0
p(15)=-14061288187707477104464845D0/35184372088832D0
p(14)=-974289645931023019776572595D0/562949953421312D0
p(13)=144483000592453610567900626155D0/4503599627370496D0
p(12)=7079170685683534011791963440605D0/18014398509481984D0
p(11)=57063078564052886214162370894557D0/36028797018963968D0
p(10)=-27650873494903018971933761124915D0/36028797018963968D0
p(9)=-641441959755400789084497279447735D0/18014398509481984D0
p(8)=-776358156363835911942964026680595D0/4503599627370496D0
p(7)=-7974397567058086827152963496557445D0/18014398509481984D0
p(6)=-11698953582630728570229643313343213D0/18014398509481984D0
p(5)=-1857674437365958001629359052983525D0/4503599627370496D0
p(4)=4809495595975287378276611244229875D0/18014398509481984D0
p(3)=26504589049384252861409184537893125D0/36028797018963968D0
p(2)=20003336218539558834627071739613125D0/36028797018963968D0
p(1)=2765140455576880316286330097421875D0/18014398509481984D0
exacteigs(1)=2**(-10)-3
exacteigs(2)=2**(-9)-3
exacteigs(3)=2**(-8)-3
exacteigs(4)=2**(-7)-3
exacteigs(5)=2**(-6)-3
exacteigs(6)=2**(-5)-3
exacteigs(7)=2**(-4)-3
exacteigs(8)=2**(-3)-3
exacteigs(9)=2**(-2)-3
exacteigs(10)=2**(-1)-3
exacteigs(11)=2**(0)-3
exacteigs(12)=2**(1)-3
exacteigs(13)=2**(2)-3
exacteigs(14)=2**(3)-3
exacteigs(15)=2**(4)-3
exacteigs(16)=2**(5)-3
exacteigs(17)=2**(6)-3
exacteigs(18)=2**(7)-3
exacteigs(19)=2**(8)-3
exacteigs(20)=2**(9)-3
WRITE(1, '(A)', advance='no')  'Test 9, '
! DSLM
ALLOCATE(er(deg),ei(deg),berr(deg))
CALL dslm(p, deg, er, ei, berr)
CALL dsort(er, ei, deg)
WRITE(1, '(ES15.2)', advance='no') dznrm2(deg,dcmplx(er,ei)-exacteigs,1)/dznrm2(deg,exacteigs,1)
WRITE(1,'(A)',advance='no') ','
DEALLOCATE(er,ei,berr)
! POLZEROS
ALLOCATE(poly(0:deg), radius(1:deg), root(1:deg), err(deg+1)) 
DO i= 0, deg
  poly(i)=dcmplx(p(i+1), 0D0)
ENDDO
CALL polzeros(deg, poly, eps, big, small, nitmax, root, radius, err, iter)
CALL zsort(root, deg)
WRITE(1, '(ES15.2)', advance='no') dznrm2(deg,root-exacteigs,1)/dznrm2(deg,exacteigs,1)
WRITE(1,'(A)', advance='no') ','
DEALLOCATE(poly,radius,root,err)
! AMVW
ALLOCATE(REIGS(deg),IEIGS(deg),ITS(deg), p2(deg))
DO i=1,deg
    p2(deg-i+1)=p(i)
ENDDO
CALL damvw(deg, p2, REIGS, IEIGS, ITS, FLAG)
CALL dsort(REIGS,IEIGS,deg)
WRITE(1, '(ES15.2)') dznrm2(deg,dcmplx(REIGS,IEIGS)-exacteigs,1)/dznrm2(deg,exacteigs,1)
DEALLOCATE(REIGS,IEIGS,ITS,p2)
! End Test 9
DEALLOCATE(exacteigs, p)

! Test 10: Chebysehv Deg 20 Polynomial
deg=20
ALLOCATE(exacteigs(deg), p(deg+1))
p(21) = 1D0
p(20) = 0D0
p(19) =-2621440D0/524288D0
p(18) = 0D0
p(17) = 5570560D0/524288D0
p(16) = 0D0
p(15) =-6553600D0/524288D0
p(14) = 0D0
p(13) = 4659200D0/524288D0
p(12) = 0D0
p(11)=-2050048D0/524288D0
p(10)= 0D0
p(9)= 549120D0/524288D0
p(8)= 0D0
p(7)=-84480D0/524288D0
p(6)= 0D0
p(5)= 6600D0/524288D0
p(4)= 0D0
p(3)=-200D0/524288D0
p(2)= 0D0
p(1)= 1D0/524288D0
exacteigs(20)  = dcmplx(cos(1D0/40D0*pi),0D0)
exacteigs(19)  = dcmplx(cos(3D0/40D0*pi),0D0)
exacteigs(18)  = dcmplx(cos(5D0/40D0*pi),0D0)
exacteigs(17)  = dcmplx(cos(7D0/40D0*pi),0D0)
exacteigs(16)  = dcmplx(cos(9D0/40D0*pi),0D0)
exacteigs(15)  = dcmplx(cos(11D0/40D0*pi),0D0)
exacteigs(14)  = dcmplx(cos(13D0/40D0*pi),0D0)
exacteigs(13)  = dcmplx(cos(15D0/40D0*pi),0D0)
exacteigs(12)  = dcmplx(cos(17D0/40D0*pi),0D0)
exacteigs(11)  = dcmplx(cos(19D0/40D0*pi),0D0)
exacteigs(10)  = dcmplx(cos(21D0/40D0*pi),0D0)
exacteigs(9)  = dcmplx(cos(23D0/40D0*pi),0D0)
exacteigs(8)  = dcmplx(cos(25D0/40D0*pi),0D0)
exacteigs(7)  = dcmplx(cos(27D0/40D0*pi),0D0)
exacteigs(6)  = dcmplx(cos(29D0/40D0*pi),0D0)
exacteigs(5)  = dcmplx(cos(31D0/40D0*pi),0D0)
exacteigs(4)  = dcmplx(cos(33D0/40D0*pi),0D0)
exacteigs(3)  = dcmplx(cos(35D0/40D0*pi),0D0)
exacteigs(2)  = dcmplx(cos(37D0/40D0*pi),0D0)
exacteigs(1)  = dcmplx(cos(39D0/40D0*pi),0D0)
WRITE(1, '(A)', advance='no')  'Test 10, '
! DSLM
ALLOCATE(er(deg),ei(deg),berr(deg))
CALL dslm(p, deg, er, ei, berr)
CALL dsort(er, ei, deg)
WRITE(1, '(ES15.2)', advance='no') dznrm2(deg,dcmplx(er,ei)-exacteigs,1)/dznrm2(deg,exacteigs,1)
WRITE(1,'(A)',advance='no') ','
DEALLOCATE(er,ei,berr)
! POLZEROS
ALLOCATE(poly(0:deg), radius(1:deg), root(1:deg), err(deg+1)) 
DO i= 0, deg
  poly(i)=dcmplx(p(i+1), 0D0)
ENDDO
CALL polzeros(deg, poly, eps, big, small, nitmax, root, radius, err, iter)
CALL zsort(root, deg)
WRITE(1, '(ES15.2)', advance='no') dznrm2(deg,root-exacteigs,1)/dznrm2(deg,exacteigs,1)
WRITE(1,'(A)', advance='no') ','
DEALLOCATE(poly,radius,root,err)
! AMVW
ALLOCATE(REIGS(deg),IEIGS(deg),ITS(deg), p2(deg))
DO i=1,deg
    p2(deg-i+1)=p(i)
ENDDO
CALL damvw(deg, p2, REIGS, IEIGS, ITS, FLAG)
CALL dsort(REIGS,IEIGS,deg)
WRITE(1, '(ES15.2)') dznrm2(deg,dcmplx(REIGS,IEIGS)-exacteigs,1)/dznrm2(deg,exacteigs,1)
DEALLOCATE(REIGS,IEIGS,ITS,p2)
! End Test 10
DEALLOCATE(exacteigs, p)

 CLOSE(1)

CALL EXECUTE_COMMAND_LINE('python testComparison.py')

CONTAINS

!********************************************************
!                       DSORT                           *
!******************************************************** 
SUBROUTINE dsort(er, ei, n)
IMPLICIT NONE
!scalar arguments
INTEGER, INTENT(IN)             :: n
!array arguments
DOUBLE PRECISION, INTENT(INOUT) :: er(*), ei(*)
!local scalars
INTEGER                         :: i, j, k
DOUBLE PRECISION                :: tr, ti

DO i=1,n
    tr=er(i); ti=ei(i); j=i
    DO k=i+1,n
        IF (tr>er(k)) THEN
            tr=er(k); ti=ei(k); j=k
        ENDIF
    ENDDO
    IF (j>i) THEN
        tr=er(i); ti=ei(i)
        er(i)=er(j); ei(i)=ei(j)
        er(j)=tr; ei(j)=ti
    ENDIF
ENDDO
END SUBROUTINE dsort

!********************************************************
!                       ZSORT                           *
!********************************************************             
SUBROUTINE zsort(root, n)
IMPLICIT NONE
!scalar arguments
INTEGER, INTENT(IN)             :: n
!array arguments
DOUBLE COMPLEX, INTENT(INOUT)   :: root(*)
!local scalars
INTEGER                         :: i, j, k
DOUBLE COMPLEX                  :: temp

DO i=1,n
    temp=root(i); j=i
    DO k=i+1,n
        IF (dble(temp)>dble(root(k))) THEN
            temp=root(k); j=k
        ENDIF
    ENDDO
    IF (j>i) THEN
        temp=root(i)
        root(i)=root(j)
        root(j)=temp
    ENDIF
ENDDO
END SUBROUTINE zsort

END PROGRAM testComparison
