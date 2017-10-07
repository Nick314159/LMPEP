PROGRAM conv_test
IMPLICIT NONE
INTEGER                                         :: deg, itmax
DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE     :: p, berr, er, ei, error
DOUBLE COMPLEX, DIMENSION(:), ALLOCATABLE       :: exacteigs
!loop variables
INTEGER                                         :: i
!intrinsic subroutines
INTRINSIC                                       :: dble, dcmplx, dimag, epsilon
!external subroutines
EXTERNAL                                        :: dslm_conv, dslm1_conv, dsam_conv
!array varaiables
CHARACTER(LEN=100)                              :: arg

itmax = 60

CALL getarg(1,arg)  
READ(arg, *) deg


!DSLM 
OPEN(UNIT=1,FILE="results/dslm_conv_test_results.txt")
ALLOCATE(berr(deg),er(deg),ei(deg),p(deg+1), exacteigs(deg), error(itmax))

!Create random problem
CALL daruv(deg+1,p)

!Solve for eigenvalues
CALL dslm(p, deg, er, ei, berr)
DO i=1,deg
  exacteigs(i) = dcmplx(er(i), ei(i))
ENDDO
DEALLOCATE(berr,er,ei)

!Test for convergence rate
ALLOCATE(berr(deg),er(deg),ei(deg))
CALL dslm_conv(p, deg, er, ei, berr, exacteigs, error)
DO i =1,itmax
  WRITE(1, '(ES15.2E3)') error(i)  
ENDDO
DEALLOCATE(berr,er,ei, exacteigs, error)

 CLOSE(UNIT=1)
!--------------------------------------

!DSLM1 
OPEN(UNIT=1,FILE="results/dslm1_conv_test_results.txt")
ALLOCATE(berr(deg),er(deg),ei(deg), exacteigs(deg), error(itmax))

!Solve for eigenvalues
CALL dslm1_conv(p, deg, er, ei, berr, exacteigs, error)
DO i=1,deg
  exacteigs(i) = dcmplx(er(i), ei(i))
ENDDO
DEALLOCATE(berr,er,ei)

!Test for convergence rate
ALLOCATE(berr(deg),er(deg),ei(deg))
CALL dslm1_conv(p, deg, er, ei, berr, exacteigs, error)
DO i =1,itmax
  WRITE(1, '(ES15.2E3)') error(i)  
ENDDO
DEALLOCATE(berr,er,ei, exacteigs, error)

 CLOSE(UNIT=1)
!--------------------------------------

!DSAM 
OPEN(UNIT=1,FILE="results/dsam_conv_test_results.txt")
ALLOCATE(berr(deg),er(deg),ei(deg), exacteigs(deg), error(itmax))

!Solve for eigenvalues
CALL dsam_conv(p, deg, er, ei, berr, exacteigs, error)
DO i=1,deg
  exacteigs(i) = dcmplx(er(i), ei(i))
ENDDO
DEALLOCATE(berr,er,ei)

!Test for convergence rate
ALLOCATE(berr(deg),er(deg),ei(deg))
CALL dsam_conv(p, deg, er, ei, berr, exacteigs, error)
DO i =1,itmax
  WRITE(1, '(ES15.2E3)') error(i)  
ENDDO
DEALLOCATE(berr,er,ei,p, exacteigs, error)

 CLOSE(UNIT=1)


CALL EXECUTE_COMMAND_LINE('python conv_test.py')
END PROGRAM conv_test
