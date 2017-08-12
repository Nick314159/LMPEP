!>\brief <b> Compute  the upper convex hull of the set </b>
!>\note Borrowed from \ref bini
!>\par Purpose:
!>\verbatim
!> Compute  the upper convex hull of the set (i,a(i)), i.e., the set of  
!> vertices (i_k,a(i_k)), k=1,2,...,m, such that the points (i,a(i)) lie 
!> below the straight lines passing through two consecutive vertices.    
!> The abscissae of the vertices of the convex hull equal the indices of 
!> the TRUE  components of the logical output vector H.                  
!> The used method requires O(nlog n) comparisons and is based on a      
!> divide-and-conquer technique. Once the upper convex hull of two       
!> contiguous sets  (say, {(1,a(1)),(2,a(2)),...,(k,a(k))} and           
!> {(k,a(k)), (k+1,a(k+1)),...,(q,a(q))}) have been computed, then       
!> the upper convex hull of their union is provided by the subroutine    
!> CMERGE. The program starts with sets made up by two consecutive       
!> points, which trivially constitute a convex hull, then obtains sets   
!> of 3,5,9... points,  up to  arrive at the entire set.                 
!> The program uses the subroutine  CMERGE; the subroutine CMERGE uses   
!> the subroutines LEFT, RIGHT and CTEST. The latter tests the convexity 
!> of the angle formed by the points (i,a(i)), (j,a(j)), (k,a(k)) in the 
!> vertex (j,a(j)) up to within a given tolerance TOLER, where i<j<k.  
!>\endverbatim
!>\param[in] n
!>\verbatim Size of the vector h \endverbatim
!>\param[in] a
!>\verbatim The set of points \endverbatim
!>\param[out] h
!>\verbatim Logical vector for which TRUE components represent the vertices of the hull of the set of points a \endverbatim
!************************************************************************
SUBROUTINE cnvex(n, a, h)
IMPLICIT NONE
INTEGER, INTENT(IN)             :: n
LOGICAL, INTENT(OUT)            :: h(*)
DOUBLE PRECISION, INTENT(IN)    :: a(*)
!intrinsic functions
INTRINSIC                       :: INT,MAX
!external subroutines
EXTERNAL                        :: cmerge
! Local variables
INTEGER                         :: i, j, k, m, nj, jc

h(1:n) = .TRUE.

! compute K such that N-2 <= 2**K < N-1
k = INT(LOG(n-2.0D0)/LOG(2.0D0))
IF(2**(k+1) <= (n-2)) k = k+1

! For each M=1,2,4,8,...,2**K, consider the NJ pairs of consecutive
! sets made up by M+1 points having the common vertex
! (JC,A(JC)), where JC=M*(2*J+1)+1 and J=0,...,NJ,
! NJ = MAX(0, INT((N-2-M)/(M+M))).
! Compute the upper convex hull of their union by means of subroutine CMERGE
m = 1
DO i = 0, k
  nj = MAX(0, INT((n-2-m)/(m+m)))
  DO j = 0, nj
    jc = (j+j+1)*m+1
    CALL cmerge(n, a, jc, m, h)
  ENDDO
  m = m+m
ENDDO
RETURN
END SUBROUTINE cnvex



