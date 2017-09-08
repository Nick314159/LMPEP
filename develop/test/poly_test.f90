PROGRAM poly_test
USE poly_zeroes
IMPLICIT NONE
INTEGER                                         :: clock, clock_rate, clock_start, clock_stop
INTEGER                                         :: i, j, nitmax, iter
REAL(KIND=dp), DIMENSION(:), ALLOCATABLE        :: radius
DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE   :: time, backward_error
REAL(KIND=dp)                                   :: eps, big, small, aux, ru, ri
COMPLEX(KIND=dp), DIMENSION(:), ALLOCATABLE     :: root, poly
LOGICAL, DIMENSION(:), ALLOCATABLE              :: err
INTEGER                                         :: it, itmax
INTEGER                                         :: deg, flag
DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE     :: berr, er, ei, p, p2, alpha
CHARACTER(LEN=100)                              :: arg
DOUBLE PRECISION, ALLOCATABLE                   :: REIGS(:), IEIGS(:), RESIDUALS(:,:)
INTEGER, ALLOCATABLE                            :: ITS(:)
DOUBLE PRECISION                                :: t2
DOUBLE COMPLEX                                  :: a, t
!intrinsic subroutines
INTRINSIC                                       :: dabs, dble, dcmplx, getarg, maxval, random_number, system_clock
INTRINSIC                                       :: epsilon, tiny, huge
!external subroutines
EXTERNAL                                        :: dslm, damvw, dseval, dzseval, drevseval, dzrevseval, daruv, init_random_seed

!Chebysehv Polynomial
deg=20
ALLOCATE(berr(deg),er(deg),ei(deg),p(deg+1))
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

CALL dslm(p, deg, er, ei, berr)
PRINT*, maxval(berr)

END PROGRAM poly_test
