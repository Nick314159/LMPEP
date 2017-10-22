!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!   Thomas R. Cameron* and Nikolas I. Steckley#
!                                                               
!   *Dept. Mathematics and Computer Science, Davidson College
!   #Dept. DiscoverOrg LLC.
!       
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!   Last modified 22 October 2017
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
!   The error vector for all polynomials is written to the
!   file: "results/testConvergence.csv".
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
PROGRAM testComparison
USE poly_zeroes
IMPLICIT NONE
! Compute variables
INTEGER :: clock, clock_rate, clock_start, clock_stop
INTEGER :: deg, i, it, itmax
! POLZERS variables
REAL(KIND=dp)   :: eps, big, small
LOGICAL, DIMENSION(:), ALLOCATABLE :: err
REAL(KIND=dp), DIMENSION(:), ALLOCATABLE :: radius
COMPLEX(KIND=dp), DIMENSION(:), ALLOCATABLE     :: root, poly


END PROGRAM testComparison
