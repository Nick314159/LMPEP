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
DOUBLE COMPLEX, DIMENSION(:), ALLOCATABLE       :: exacteigs
CHARACTER(LEN=100)                              :: arg
DOUBLE PRECISION, ALLOCATABLE                   :: REIGS(:), IEIGS(:), RESIDUALS(:,:)
INTEGER, ALLOCATABLE                            :: ITS(:)
DOUBLE PRECISION                                :: t2
DOUBLE COMPLEX                                  :: a, t
DOUBLE PRECISION, PARAMETER                     :: pi = 3.14159265358979323d0
!intrinsic subroutines
INTRINSIC                                       :: dabs, dble, dcmplx, getarg, maxval, random_number, system_clock
INTRINSIC                                       :: epsilon, tiny, huge
!external subroutines
EXTERNAL                                        :: dslm, damvw, dseval, dzseval, drevseval, dzrevseval, daruv, init_random_seed
!external functions
DOUBLE PRECISION                                :: dzmod
EXTERNAL                                        :: dzmod

!Chebysehv Polynomial
deg=20
ALLOCATE(exacteigs(deg),p(deg+1))
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
exacteigs(20)  = complex(cos(1d0/40d0*pi),0d0)
exacteigs(19)  = complex(cos(3d0/40d0*pi),0d0)
exacteigs(18)  = complex(cos(5d0/40d0*pi),0d0)
exacteigs(17)  = complex(cos(7d0/40d0*pi),0d0)
exacteigs(16)  = complex(cos(9d0/40d0*pi),0d0)
exacteigs(15)  = complex(cos(11d0/40d0*pi),0d0)
exacteigs(14)  = complex(cos(13d0/40d0*pi),0d0)
exacteigs(13)  = complex(cos(15d0/40d0*pi),0d0)
exacteigs(12)  = complex(cos(17d0/40d0*pi),0d0)
exacteigs(11)  = complex(cos(19d0/40d0*pi),0d0)
exacteigs(10)  = complex(cos(21d0/40d0*pi),0d0)
exacteigs(9)  = complex(cos(23d0/40d0*pi),0d0)
exacteigs(8)  = complex(cos(25d0/40d0*pi),0d0)
exacteigs(7)  = complex(cos(27d0/40d0*pi),0d0)
exacteigs(6)  = complex(cos(29d0/40d0*pi),0d0)
exacteigs(5)  = complex(cos(31d0/40d0*pi),0d0)
exacteigs(4)  = complex(cos(33d0/40d0*pi),0d0)
exacteigs(3)  = complex(cos(35d0/40d0*pi),0d0)
exacteigs(2)  = complex(cos(37d0/40d0*pi),0d0)
exacteigs(1)  = complex(cos(39d0/40d0*pi),0d0)


!LMPEP
ALLOCATE(berr(deg),er(deg),ei(deg))

CALL dslm(p, deg, er, ei, berr)
PRINT*, maxval(berr)

CALL dsort(er, ei, deg)
PRINT*, 'LMPEP Absolute Error: Deg 20 Chebyshev Poly'
DO i=1,deg
    PRINT*, dzmod(dble(exacteigs(i))-er(i),dimag(exacteigs(i))-ei(i))
ENDDO
PRINT*, ''

!POLZEROS
ALLOCATE(poly(0:deg), radius(1:deg), root(1:deg), err(deg+1))

eps    = epsilon(1.0D0)
small  = tiny(1.0D0)
big    = huge(1.0D0)

DO i= 0, deg
    poly(i)=dcmplx(p(i+1), 0.0D0)
END DO

CALL polzeros(deg, poly, eps, big, small, nitmax, root, radius, err, iter)

CALL zsort(root, deg)
PRINT*, 'POLZEROS Absolute Error: Deg 20 Chebyshev Poly'
DO i=1,deg
    PRINT*, dzmod(dble(exacteigs(i))-dble(root(i)),dimag(exacteigs(i))-dimag(root(i)))
ENDDO
PRINT*, ''

CONTAINS

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

END PROGRAM poly_test
