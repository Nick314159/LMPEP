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
INTEGER :: deg, i, itmax, nmax
DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE :: p
DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE :: error

! Allocate error vector
itmax = 60 ! number of iterations allowed in dslm
nmax = 2   ! number of tests
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
WRITE(1, '(A)') "Deg 20 Chebyshev Poly"
CALL test(deg, p, error(:,1))
DEALLOCATE(p)

DO i=1,10
    WRITE(*,*) error(i,1)
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
