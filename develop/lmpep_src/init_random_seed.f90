!>\author Thomas R. Cameron, Davidson College
!>\author Nikolas I. Steckley, Portland State University
!>\date 2017
!>\brief <b> Init Random Seed initialzes the random number generator seed </b>
!>\par Purpose:
!>\verbatim
!> Init Random Seed sets the random number gernerater's seed, ensuring truely random problems each run of the program
!>\endverbatim
!***********************************************************************
SUBROUTINE init_random_seed()
INTEGER                             :: i, n, clock
INTEGER, DIMENSION(:), ALLOCATABLE  :: seed
!intrinsic procedures
INTRINSIC                           :: count, random_seed, size, system_clock

CALL random_seed(size = n)
ALLOCATE(seed(n))

CALL system_clock(count=clock)
seed = clock + 37 * (/ (i - 1, i = 1, n) /)
CALL random_seed(PUT = seed)

DEALLOCATE(seed)
RETURN
END SUBROUTINE init_random_seed
