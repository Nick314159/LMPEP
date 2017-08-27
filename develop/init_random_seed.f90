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

CALL RANDOM_SEED(size = n)
ALLOCATE(seed(n))

CALL SYSTEM_CLOCK(COUNT=clock)
seed = clock + 37 * (/ (i - 1, i = 1, n) /)
CALL RANDOM_SEED(PUT = seed)

DEALLOCATE(seed)
END SUBROUTINE init_random_seed
