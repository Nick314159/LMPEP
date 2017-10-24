!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!   Thomas R. Cameron* and Nikolas I. Steckley#
!                                                               
!   *Dept. Mathematics and Computer Science, Davidson College
!   #Dept. DiscoverOrg LLC.
!       
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!   Last modified 22 October 2017
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!   Test Convergence on several famous polynomials 
!   
!   The root approximations are computed using DSLM and then
!   an error vector is computed using DSLM_CONV.
!   
!   The error vector for all polynomials is written to the
!   file: "results/testConvergence.csv".
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
PROGRAM testConvergence
IMPLICIT NONE
! Compute Variables
INTEGER :: deg, i, j, itmax, ntest
DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE :: p
DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE :: error

! Allocate error vector
itmax = 60 ! number of iterations allowed in dslm
ntest = 12   ! number of tests
ALLOCATE(error(itmax+1,ntest))

! Open file
OPEN(UNIT=1,FILE="results/testConvergence.csv")

! Test 1: Wilkinson Deg 10 Polynomial
deg=10
ALLOCATE(p(deg+1))
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
WRITE(1, '(A)', advance='no')  'Test 1'
WRITE(1,'(A)', advance='no') ','
CALL test(deg, p, error(:,1))
DEALLOCATE(p)

! Test 2: Wilkinson Deg 15 Polynomial
deg=15
ALLOCATE(p(deg+1))
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
WRITE(1, '(A)', advance='no')  'Test 2'
WRITE(1,'(A)', advance='no') ','
CALL test(deg, p, error(:,2))
DEALLOCATE(p)

! Test 3: Wilkinson Deg 20 Polynomial
deg=20
ALLOCATE(p(deg+1))
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
WRITE(1, '(A)', advance='no')  'Test 3'
WRITE(1,'(A)', advance='no') ','
CALL test(deg, p, error(:,3))
DEALLOCATE(p)

! Test 4: scaled and shifted Wilkinson Deg 20 Polynomial
deg=20
ALLOCATE(p(deg+1))
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
WRITE(1, '(A)', advance='no')  'Test 4'
WRITE(1,'(A)', advance='no') ','
CALL test(deg, p, error(:,4))
DEALLOCATE(p)


! Test 5: reverse Wilkinson Deg 10 polynomial
deg=10
ALLOCATE(p(deg+1))
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
WRITE(1, '(A)', advance='no')  'Test 5'
WRITE(1,'(A)', advance='no') ','
CALL test(deg, p, error(:,5))
DEALLOCATE(p)

! Test 6: reverse Wilkinson Deg 15 polynomial
deg=15
ALLOCATE(p(deg+1))
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
WRITE(1, '(A)', advance='no')  'Test 6'
WRITE(1,'(A)', advance='no') ','
CALL test(deg, p, error(:,6))
DEALLOCATE(p)

! Test 7: reverse Wilkinson Deg 20 polynomial
deg=20
ALLOCATE(p(deg+1))
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
WRITE(1, '(A)', advance='no')  'Test 7'
WRITE(1,'(A)', advance='no') ','
CALL test(deg, p, error(:,7))
DEALLOCATE(p)

! Test 8: prescribed roots of varying scale
deg=20
ALLOCATE(p(deg+1))
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
WRITE(1, '(A)', advance='no')  'Test 8'
WRITE(1,'(A)', advance='no') ','
CALL test(deg, p, error(:,8))
DEALLOCATE(p)

! Test 9: prescribed roots of varying scale - 3
deg=20
ALLOCATE(p(deg+1))
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
WRITE(1, '(A)', advance='no')  'Test 9'
WRITE(1,'(A)', advance='no') ','
CALL test(deg, p, error(:,9))
DEALLOCATE(p)

! Test 10: Chebysehv Deg 20 Polynomial
deg=20
ALLOCATE(p(deg+1))
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
WRITE(1, '(A)', advance='no')  'Test 10'
WRITE(1,'(A)', advance='no') ','
CALL test(deg, p, error(:,10))
DEALLOCATE(p)

! Test 11: Deg 20 polynomial whose coefficients are all 1
deg=20
ALLOCATE(p(deg+1))
p(1:deg+1)=1D0
WRITE(1, '(A)', advance='no')  'Test 11'
WRITE(1, '(A)', advance='no') ','
CALL test(deg, p, error(:,11))
DEALLOCATE(p)

! Test 12: Deg 20 Random Polynomial
deg=20
ALLOCATE(p(deg+1))
CALL daruv(deg+1,p)
WRITE(1, '(A)') 'Random Poly'
CALL test(deg, p, error(:,12))
DEALLOCATE(p)

DO i=1,10
  DO j=1,ntest
    WRITE(1,'(ES15.2)', advance='no') error(i,j)
    IF (j.NE.ntest) THEN
      WRITE(1,'(A)', advance='no') ','
    ELSE
      WRITE(1,*)
    ENDIF
  ENDDO
ENDDO

! Deallocate error vector
DEALLOCATE(error)

! Close results file
 CLOSE(1)

CALL EXECUTE_COMMAND_LINE('python testConvergence.py')


CONTAINS

!********************************************************
!                       TEST                            *
!******************************************************** 
SUBROUTINE test(deg,p,error)
IMPLICIT NONE
! scalar arguments
INTEGER, INTENT(IN) :: deg
! array arguments
DOUBLE PRECISION, INTENT(IN) :: p(*)
DOUBLE PRECISION, INTENT(INOUT) :: error(*)
! DSLM variables
DOUBLE PRECISION, DIMENSION(deg) :: berr, er, ei, rr, ri


! Call dslm and store results
CALL dslm(p, deg, er, ei, berr)
DO i=1,deg
    rr(i)=er(i)
    ri(i)=ei(i)
ENDDO
! Call dslm_conv
CALL dslm_conv(p, deg, er, ei, berr, rr, ri, error)

END SUBROUTINE


END PROGRAM testConvergence
