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

!Wilkonson Polynomial
p(21)= 243290200817664D4
p(20)= -875294803676160D4
p(19)= 138037597536407D5
p(18)= -128709312451510D5
p(17)= 803781182264505D4
p(16)= -359997951794761D4
p(15)= 120664780378037D4
p(14)= -311333643161391D3
p(13)= 630308120992949D2
p(12)= -101422998655115D2
p(11)= 130753501054049D2
p(10)= -135585182899530D0
p(9) = 11310276995381D0
p(8) = -756111184500D0
p(7) = 40171771630D0
p(6) = -1672280820D0
p(5) = 53327946D0
p(4) = -1256850D0
p(3) = 20615D0
p(2) = -210D0
p(1) = 1D0
exacteigs(20)  = 20
exacteigs(19)  = 19
exacteigs(18)  = 18
exacteigs(17)  = 17
exacteigs(16)  = 16
exacteigs(15)  = 15
exacteigs(14)  = 14
exacteigs(13)  = 13
exacteigs(12)  = 12
exacteigs(11)  = 11
exacteigs(10)  = 10
exacteigs(9)  = 9
exacteigs(8)  = 8
exacteigs(7)  = 7
exacteigs(6)  = 6
exacteigs(5)  = 5
exacteigs(4)  = 4
exacteigs(3)  = 3
exacteigs(2)  = 2
exacteigs(1)  = 1




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
DEALLOCATE(poly, radius, root, err)

!AMVW
ALLOCATE(REIGS(deg),IEIGS(deg),ITS(deg), p2(deg))
DO i=1,deg
    p2(deg-i+1)=p(i)/p(deg+1)
ENDDO

CALL damvw(deg, p2, REIGS, IEIGS, ITS, FLAG)
CALL dsort(REIGS, IEIGS, deg)
PRINT*, 'AMVW Absolute Error: Deg 20 Chebyshev Poly'
DO i=1,deg
    PRINT*, dzmod(dble(exacteigs(i))-REIGS(i),dimag(exacteigs(i))-IEIGS(i))
ENDDO
DEALLOCATE(p, p2, REIGS,IEIGS,ITS, berr)
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
