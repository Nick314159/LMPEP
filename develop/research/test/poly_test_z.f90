PROGRAM poly_test
IMPLICIT NONE
INTEGER                                         :: deg
DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE     :: p
DOUBLE COMPLEX, DIMENSION(:), ALLOCATABLE       :: exacteigs
DOUBLE PRECISION, PARAMETER                     :: pi = 3.1415926535897932d0
!intrinsic subroutines
INTRINSIC                                       :: cos, dcmplx

OPEN(UNIT=1,FILE="results/poly_test_z_results.txt")

!Cylcotomic Deg 20 Polynomial
deg=20
ALLOCATE(exacteigs(deg), p(deg+1))
p(1)=1D0
p(2)=0D0
p(3)=0D0
p(4)=0D0
p(5)=0D0
p(6)=1D0
p(7)=0D0
p(8)=0D0
p(9)=0D0
p(10)=0D0
p(11)=1D0
p(12)=0D0
p(13)=0D0
p(14)=0D0
p(15)=0D0
p(16)=1D0
p(17)=0D0
p(18)=0D0
p(19)=0D0
p(20)=0D0
p(21)=1D0
exacteigs(1)=dcmplx(-0.992114701314478, -0.125333233564304)
exacteigs(2)=dcmplx(0.968583161128631, 0.248689887164855)
exacteigs(3)=dcmplx(-0.929776485888251, -0.368124552684678)
exacteigs(4)=dcmplx(0.876306680043864, 0.481753674101715)
exacteigs(5)=dcmplx(0.728968627421412, 0.684547105928689)
exacteigs(6)=dcmplx(-0.637423989748690, -0.770513242775789)
exacteigs(7)=dcmplx(0.535826794978997, 0.844327925502015)
exacteigs(8)=dcmplx(-0.425779291565073, -0.904827052466020)
exacteigs(9)=dcmplx(-0.187381314585725, -0.982287250728689)
exacteigs(10)=dcmplx(0.062790519529313, 0.998026728428272)
exacteigs(11)=dcmplx(0.062790519529313, -0.998026728428272)
exacteigs(12)=dcmplx(-0.187381314585725, 0.982287250728689)
exacteigs(13)=dcmplx(-0.425779291565073, 0.904827052466020)
exacteigs(14)=dcmplx(0.535826794978997, -0.844327925502015)
exacteigs(15)=dcmplx(-0.637423989748690, 0.770513242775789)
exacteigs(16)=dcmplx(0.728968627421412, -0.684547105928689)
exacteigs(17)=dcmplx(0.876306680043864, -0.481753674101715)
exacteigs(18)=dcmplx(-0.929776485888251, 0.368124552684678)
exacteigs(19)=dcmplx(0.968583161128631, -0.248689887164855)
exacteigs(20)=dcmplx(-0.992114701314478, 0.125333233564304)
WRITE(1, '(A)')  "Deg 20 Cyclotomic Poly"
CALL test(deg, p, exacteigs)
DEALLOCATE(exacteigs, p)

 CLOSE(UNIT=1)
CONTAINS

!********************************************************
!                       TEST                            *
!******************************************************** 
SUBROUTINE test(deg, p, exacteigs)
USE poly_zeroes
IMPLICIT NONE
!scalar arguments
INTEGER, INTENT(IN)                             :: deg
!array arguments
DOUBLE PRECISION, INTENT(IN)                    :: p(*) 
DOUBLE COMPLEX, INTENT(IN)                      :: exacteigs(*)
!LMPEP variables
DOUBLE PRECISION, DIMENSION(deg)                :: berr, er, ei
!POLZEROS variables
INTEGER                                         :: iter, nitmax
LOGICAL, DIMENSION(deg+1)                       :: err
REAL(KIND=dp), DIMENSION(1:deg)                 :: radius
COMPLEX(KIND=dp), DIMENSION(0:deg)              :: poly
COMPLEX(KIND=dp), DIMENSION(1:deg)              :: root
REAL(KIND=dp), PARAMETER                        :: eps=epsilon(1.0D0), big=huge(1.0D0), small=tiny(1.0D0)
!AMVW variables
INTEGER                                         :: FLAG
INTEGER, DIMENSION(deg)                         :: ITS
DOUBLE PRECISION, DIMENSION(deg)                :: REIGS, IEIGS, P2
!loop variables
INTEGER                                         :: i
!intrinsic subroutines
INTRINSIC                                       :: dble, dcmplx, dimag, epsilon, tiny, huge
!external subroutines
EXTERNAL                                        :: dslm, damvw
!external functions
DOUBLE PRECISION                                :: dzmod
EXTERNAL                                        :: dzmod


!LMPEP
CALL dslm(p, deg, er, ei, berr)
CALL dsort(er, ei, deg)
WRITE(1, '(A)') 'LMPEP Absolute Error:'
DO i=1,deg
   WRITE(1, *) dzmod(dble(exacteigs(i))-er(i),dimag(exacteigs(i))-ei(i))
ENDDO
PRINT*, ''

!POLZEROS
nitmax=60
DO i=0,deg
    poly(i)=dcmplx(p(i+1), 0.0D0)
ENDDO
CALL polzeros(deg, poly, eps, big, small, nitmax, root, radius, err, iter)
CALL zsort(root, deg)
WRITE(1, '(A)') 'POLZEROS Absolute Error:'
DO i=1,deg
    WRITE(1, *) dzmod(dble(exacteigs(i))-dble(root(i)),dimag(exacteigs(i))-dimag(root(i)))
ENDDO
PRINT*, ''

!AMVW
FLAG=1
DO i=1,deg
    p2(deg-i+1)=p(i)/p(deg+1)
ENDDO
CALL damvw(deg, P2, REIGS, IEIGS, ITS, FLAG)
CALL dsort(REIGS, IEIGS, deg)
WRITE(1, '(A)') 'AMVW Absolute Error:'
DO i=1,deg
    WRITE(1, *) dzmod(dble(exacteigs(i))-REIGS(i),dimag(exacteigs(i))-IEIGS(i))
ENDDO
PRINT*, ''

END SUBROUTINE test

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
!SUBROUTINE zsort(root, n)
!IMPLICIT NONE
!scalar arguments
!INTEGER, INTENT(IN)             :: n
!array arguments
!DOUBLE COMPLEX, INTENT(INOUT)   :: root(*)
!local scalars
!INTEGER                         :: i, j, k
!DOUBLE COMPLEX                  :: temp

!DO i=1,n
!    temp=root(i); j=i
!    DO k=i+1,n
!        IF (dble(temp)>dble(root(k))) THEN
!            temp=root(k); j=k
!        ENDIF
!    ENDDO
!    IF (j>i) THEN
!        temp=root(i)
!        root(i)=root(j)
!        root(j)=temp
!    ENDIF
!ENDDO
!END SUBROUTINE zsort

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
        IF (atan2(aimag(temp), real(temp))>atan2(aimag(root(k)), real(root(k)))) THEN
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
