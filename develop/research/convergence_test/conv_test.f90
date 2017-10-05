PROGRAM conv_test
IMPLICIT NONE
INTEGER, PARAMETER     :: dp = SELECTED_REAL_KIND(15, 60)
INTEGER                                         :: deg
DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE     :: p, berr, er, ei
DOUBLE COMPLEX, DIMENSION(:), ALLOCATABLE       :: exacteigs
DOUBLE PRECISION, PARAMETER                     :: pi = 3.1415926535897932d0
!loop variables
INTEGER                                         :: i, j
!intrinsic subroutines
INTRINSIC                                       :: dble, dcmplx, dimag, epsilon, tiny, huge
!external subroutines
EXTERNAL                                        :: dslm_conv
!external functions
DOUBLE PRECISION                                :: dzmod
EXTERNAL                                        :: dzmod
!scalar variables
INTEGER                                         :: itmax
!array varaiables
DOUBLE COMPLEX, DIMENSION(:), ALLOCATABLE       :: error
CHARACTER(LEN=100)                              :: arg

CALL getarg(1,arg)  
READ(arg, *) deg
CALL getarg(2,arg)
READ(arg, *) itmax

OPEN(UNIT=1,FILE="results/conv_test_results.txt")

!DSLM 
!Create random problem
ALLOCATE(berr(deg),er(deg),ei(deg),p(deg+1), exacteigs(deg), error(itmax))
CALL daruv(deg+1,p)

!Solve for eigenvalues
CALL dslm(p, deg, er, ei, berr)
DO i=1,deg
  exacteigs(i) = dcmplx(er(i), ei(i))
ENDDO
DEALLOCATE(berr,er,ei)

!Test for convergence rate
ALLOCATE(berr(deg),er(deg),ei(deg))
CALL dslm_conv(p, deg, er, ei, berr, exacteigs, error, itmax)
CALL dsort(er, ei, deg)
DO i =1,itmax
  WRITE(1, '(ES15.2)') error(i)  
ENDDO
DEALLOCATE(berr,er,ei,p, exacteigs, error)

 CLOSE(UNIT=1)
CALL EXECUTE_COMMAND_LINE('python conv_test.py')
CONTAINS

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

END PROGRAM conv_test
