!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!   Thomas R. Cameron* and Nikolas I. Steckley#
!                                                               
!   *Dept. Mathematics and Computer Science, Davidson College
!   #Dept. DiscoverOrg LLC.
!       
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!   Last modified 21 October 2017
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!   Test Convergence on several famous polynomials and one random
!   polynomial of degree 100. 
!   
!   The root approximations are computed using DSLM and then
!   an error vector is computed using DSLM_CONV.
!   
!   The error vector for all polynomials is written to a the
!   file: "".
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
PROGRAM testConvergence
IMPLICIT NONE
! Compute Variables
INTEGER :: deg, i, j, itmax, nmax
DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE :: p
DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE :: error

! Allocate error vector
itmax = 60 ! number of iterations allowed in dslm
nmax = 3   ! number of tests
ALLOCATE(error(itmax+1,nmax))

! Open file
OPEN(UNIT=1,FILE="results/testConvergence.csv")

! Test 1: Chebysehv Deg 20 Polynomial
deg=20
ALLOCATE(p(deg+1))
p(21) = 1d0
p(20) = 0d0
p(19) =-2621440d0/524288d0
p(18) = 0d0
p(17) =+5570560d0/524288d0
p(16) = 0d0
p(15) =-6553600d0/524288d0
p(14) = 0d0
p(13) =+4659200d0/524288d0
p(12) = 0d0
p(11)=-2050048d0/524288d0
p(10)= 0d0
p(9)=+549120d0/524288d0
p(8)= 0d0
p(7)=-84480d0/524288d0
p(6)= 0d0
p(5)=+6600d0/524288d0
p(4)= 0d0
p(3)=-200d0/524288d0
p(2)= 0d0
p(1)= 1d0/524288d0
WRITE(1, '(A)', advance='no') "Deg 20 Chebyshev"
WRITE(1,'(A)', advance='no') ','
CALL test(deg, p, error(:,1))
DEALLOCATE(p)


!Test 2: Wilkinson Deg 20 Polynomial
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
WRITE(1, '(A)', advance='no')  "Deg 20 Wilkins"
WRITE(1,'(A)', advance='no') ','
CALL test(deg, p, error(:,2))
DEALLOCATE(p)


!Test 3: Create random polynomial
deg=20
ALLOCATE(p(deg+1))
CALL daruv(deg+1,p)
WRITE(1, '(A)') "Deg 20 Random"
CALL test(deg, p, error(:,3))
DEALLOCATE(p)

DO i=1,10
  DO j=1,nmax
    WRITE(1,'(ES15.2)', advance='no') error(i,j)
    IF (j.NE.nmax) THEN
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
