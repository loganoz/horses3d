!
!////////////////////////////////////////////////////////////////////////
!
!   @File: InterpolationAndDerivatives.f90
!   @Author: David Kopriva 
!   @Created: 2009-12-15 15:36:24 -0500 
!   @Last revision date: Mon Sep  6 22:45:02 2021
!   @Last revision author: Wojciech Laskowski (wj.laskowski@upm.es)
!   @Last revision commit: 3334a040b8cdf3201850a2deec9950c84f2dc21f
!
!//////////////////////////////////////////////////////
!
!      Contains:
!               
!               SUBROUTINE InterpolatingPolynomialVector( x, N, nodes, weights, p )
!               REAL(KIND=RP) FUNCTION EvaluateLagrangePolyDerivative( j, x, N, nodes)
!
!               ALGORITHM 30: SUBROUTINE BarycentricWeights( N, x, w )
!               ALGORITHM 31: REAL(KIND=RP) FUNCTION LagrangeInterpolation( x, N, nodes, values, weights)
!               ALGORITHM 32: SUBROUTINE PolynomialInterpolationMatrix( N, M, oldNodes, weights, newNodes, T)
!               ALGORITHM 33: SUBROUTINE InterpolateToNewPoints( N, M, T, f, fInterp )
!               ALGORITHM 34: REAL(KIND=RP) FUNCTION LagrangeInterpolatingPolynomial( j, x, N, nodes )
!               ALGORITHM 35: 
!               ALGORITHM 36: REAL(KIND=RP) FUNCTION LagrangeInterpolantDerivative( x, N, nodes, values, weights)
!               ALGORITHM 37: SUBROUTINE PolynomialDerivativeMatrix( N, nodes, D )
!               ALGORITHM 38: SUBROUTINE mthPolynomialDerivativeMatrix( m, N, nodes, D )
!               ALGORITHM 39: SUBROUTINE EOMatrixDerivative( u, uDeriv, D, N )
!
!               The following are combined:
!                  ALGORITHM 19: MXVDerivative
!                  ALGORITHM 106: TransposeMatrixMultiply
!               as SUBROUTINE MatrixMultiplyDeriv( f, fDeriv, D, N, transp )
!
!      
!////////////////////////////////////////////////////////////////////////
!
!
!  ******
   MODULE PolynomialInterpAndDerivsModule
!  ******
!
      USE SMConstants
      IMPLICIT NONE
      
      private  
      
      public   MXV_DIRECT, MXV_TRANSPOSE
      public   Interpolator_t
      
      public   EvaluateLagrangePolyDerivative
      public   Interpolate3D, Create3DInterpolationMatrix
      public   Create3DRestrictionMatrix
      public   BarycentricWeights, PolynomialDerivativeMatrix
      public   InterpolatingPolynomialVector, PolyDerivativeVector
      public   PolynomialInterpolationMatrix
      public   MatrixMultiplyDeriv
      public   ComputeVandermonde, JacobiPolynomial, ComputeModalForm

      INTEGER, PARAMETER :: MXV_DIRECT = 1, MXV_TRANSPOSE = 2
      
      ! Interpolator type that is used in multigrid and in the plotter
      TYPE Interpolator_t
         LOGICAL                     :: Created = .FALSE.
         REAL(KIND=RP), ALLOCATABLE  :: Mat(:,:)
      END TYPE Interpolator_t
!
!    ========
     CONTAINS
!    ========
!
!    /////////////////////////////////////////////////////////////////
!
      REAL(KIND=RP) FUNCTION LagrangeInterpolatingPolynomial( j, x, N, nodes ) RESULT(p)
!
!---------------------------------------------------------------------
! Compute L_j(x) of degree N whose zeros are at the nodes using direct
! evaluation.
!---------------------------------------------------------------------
!
      INTEGER                      , INTENT(IN) :: j     !! Which polynomial
      INTEGER                      , INTENT(IN) :: N     !! Polynomial order
      REAL(KIND=RP)                , INTENT(IN) :: x     !! evaluation point
      REAL(KIND=RP), DIMENSION(0:N), INTENT(IN) :: nodes !! Interpolation nodes
!
!     ---------------
!     Local Variables
!     ---------------
!
      INTEGER :: k
      
      IF ( j == 0 )     THEN
         p = (x - nodes(1))/(nodes(j) - nodes(1))
         DO k = 2, N 
            p = p*(x - nodes(k))/(nodes(j) - nodes(k))
         END DO
      ELSE
         p = (x - nodes(0))/(nodes(j) - nodes(0))
         DO k = 1, j-1 
            p = p*(x - nodes(k))/(nodes(j) - nodes(k))
         END DO
         DO k = j+1, N
            p = p*(x - nodes(k))/(nodes(j) - nodes(k))
         END DO
      END IF
   END FUNCTION LagrangeInterpolatingPolynomial
!
!    /////////////////////////////////////////////////////////////////
!
      SUBROUTINE InterpolatingPolynomialVector( x, N, nodes, weights, p )
!
!---------------------------------------------------------------------
! Compute L_j(x), j = 0, ..., N of degree N whose zeros are at the nodes 
! using barycentric form.
!---------------------------------------------------------------------
!
      use Utilities, only: almostEqual
      INTEGER                      , INTENT(IN)  :: N       !! Polynomial order
      REAL(KIND=RP)                , INTENT(IN)  :: x       !! evaluation point
      REAL(KIND=RP), DIMENSION(0:N), INTENT(IN)  :: nodes   !! Interpolation nodes
      REAL(KIND=RP), DIMENSION(0:N), INTENT(IN)  :: weights !! Barycentric weights
      REAL(KIND=RP), DIMENSION(0:N), INTENT(INOUT) :: p
!
!     ---------------
!     Local Variables
!     ---------------
!
      INTEGER       :: j
      LOGICAL       :: xMatchesNode
      REAL(KIND=RP) :: d, t
!
!     -------------------------------------
!     See if the evaluation point is a node
!     -------------------------------------
!
      xMatchesNode = .false.
      DO j = 0, N 
         p(j) = 0.0_RP
         IF( AlmostEqual( x, nodes(j) ) )     THEN
            p(j) = 1.0_RP
            xMatchesNode = .true.
         END IF
      END DO
      IF( xMatchesNode )     RETURN
!
!     ------------------------------
!     Evaluation point is not a node
!     ------------------------------
!
      d = 0.0_RP
      DO j = 0, N 
         t = weights(j)/( x - nodes(j) )
         p(j) = t
         d = d + t
      END DO
      DO j = 0, N 
         p(j) = p(j)/d
      END DO
      
   END SUBROUTINE InterpolatingPolynomialVector
!
! /////////////////////////////////////////////////////////////////////
!
!  ---------------------------------------------------------
!  Algorithm to obtain the polynomial derivative vector in x
!  -> TODO: See if there's a more efficient way of computing
!           this with the barycentric formula.
!  ---------------------------------------------------------
   subroutine PolyDerivativeVector( x, N, nodes, p )
      implicit none
      !----------------------------------------------
      real(kind=RP)                , intent(in)    :: x     !<  Where to evaluate the derivative
      integer                      , intent(in)    :: N     !<  Polynomial order
      real(kind=RP), dimension(0:N), intent(in)    :: nodes !<  Node position
      real(kind=RP), dimension(0:N), intent(inout) :: p     !<> Poly derivative vector
      !----------------------------------------------
      integer :: j
      !----------------------------------------------
      
      do j = 0, N
         p(j) = EvaluateLagrangePolyDerivative( j, x, N, nodes)
      end do
   end subroutine PolyDerivativeVector
!
! /////////////////////////////////////////////////////////////////////
!
!---------------------------------------------------------------------
!!    Compute at x the derivative of the lagrange polynomial L_j(x) of
!!    degree N whose zeros are given by the nodes(i)
!---------------------------------------------------------------------
!
      REAL(KIND=RP) FUNCTION EvaluateLagrangePolyDerivative( j, x, N, nodes)
!
!     ---------
!     Arguments
!     ---------
!
      INTEGER, INTENT(IN)                       :: j !! Which poly
      INTEGER, INTENT(IN)                       :: N !! Order
      REAL(KIND=RP), DIMENSION(0:N), INTENT(IN) :: nodes !! Nodes
      REAL(KIND=RP)                , INTENT(IN) :: x !! Eval Pt
!
!     ---------------
!     Local Variables
!     ---------------
!
      INTEGER       :: l, m
      REAL(KIND=RP) :: hp, poly
!                                                                       
      hp = 0.0_RP
      DO l = 0,N 
         IF(l == j)     CYCLE
         poly = 1.0_RP
         DO m = 0,N 
            IF (m == l)     CYCLE
            IF (m == j)     CYCLE 
            poly = poly*(x - nodes(m))/(nodes(j) - nodes(m))
         END DO
         hp = hp + poly/(nodes(j) - nodes(l)) 
      END DO
      EvaluateLagrangePolyDerivative = hp
!                                                                       
      END FUNCTION EvaluateLagrangePolyDerivative
!
!////////////////////////////////////////////////////////////////////////
!
!     ----------------------------------------------------------------
!!    Compute the barycentric weights for polynomial interpolation
!     ----------------------------------------------------------------
!
      SUBROUTINE BarycentricWeights( N, x, w )
!
!     ---------
!     Arguments
!     ---------
!
      INTEGER :: N
      REAL(KIND=RP), DIMENSION(0:N), INTENT(IN)  :: x
      REAL(KIND=RP), DIMENSION(0:N), INTENT(OUT) :: w
!
!     ---------------
!     Local Variables
!     ---------------
!
      INTEGER :: j, k
!      
      w = 1.0_RP
      DO j = 1, N
         DO k = 0, j-1
            w(k) = w(k)*(x(k) - x(j))
            w(j) = w(j)*(x(j) - x(k))
         END DO
      END DO
      w = 1.0_RP/w

      END SUBROUTINE BarycentricWeights
!
!     ////////////////////////////////////////////////////////////////
!
!     ------------------------------------------------------------------
!!    Compute the value of the interpolant using the barycentric formula
!     ------------------------------------------------------------------
!
      REAL(KIND=RP) FUNCTION LagrangeInterpolation( x, N, nodes, values, weights)
!
!     ---------
!     Arguments
!     ---------
!
      use Utilities, only: almostEqual
      REAL(KIND=RP)                 :: x
      INTEGER                       :: N
      REAL(KIND=RP), DIMENSION(0:N) :: nodes
      REAL(KIND=RP), DIMENSION(0:N) :: values
      REAL(KIND=RP), DIMENSION(0:N) :: weights
!
!     ---------------
!     Local Variables
!     ---------------
!
      INTEGER       :: j
      REAL(KIND=RP) :: t, numerator, denominator
      
      numerator   = 0.0_RP
      denominator = 0.0_RP
      DO j = 0, N
         IF( AlmostEqual( x, nodes(j) ) )    THEN
            LagrangeInterpolation = values(j)
            RETURN 
         END IF 
         t = weights(j)/( x - nodes(j) )
         numerator = numerator + t*values(j)
         denominator = denominator + t
      END DO
      LagrangeInterpolation = numerator/denominator

      END FUNCTION LagrangeInterpolation
!
!     ////////////////////////////////////////////////////////////////
!
!     --------------------------------------------------------------------------
!!    Compute the value of the interpolant WITHOUT using the barycentric formula
!     --------------------------------------------------------------------------
!
   FUNCTION LagrangeInterpolationNoBar( x, N, nodes, j) RESULT(l)
!
!     ---------
!     Arguments
!     ---------
!
      use Utilities, only: almostEqual
      REAL(KIND=RP)                 :: l     !>  Lagrange interpolant
      REAL(KIND=RP)                 :: x     !<  Point of evaluation of interpolant
      INTEGER                       :: N     !<  Polynomial order  
      REAL(KIND=RP), DIMENSION(0:N) :: nodes !<  Nodes of Lagrange interpolation
      INTEGER                       :: j     !<  Index of polynomial to be found
!
!     ---------------
!     Local Variables
!     ---------------
!
      INTEGER                       :: i
      REAL(KIND=RP)                 :: numerator, denominator
      REAL(KIND=RP), DIMENSION(0:N) :: values
      !-----------------------------------------------------------------------------
      
      values      = 0.0_RP
      values(j)   = 1.0_RP
      numerator   = 1.0_RP
      denominator = 1.0_RP

      DO i = 0, N
         IF( AlmostEqual( x, nodes(i) ) )    THEN
            l = values(i)
            RETURN 
         ELSE IF (j.ne.i) THEN
         numerator   = numerator*(x - nodes(i))    
         denominator = denominator*(nodes(j) - nodes(i))
         END IF 
      END DO
      l = numerator/denominator

   END FUNCTION LagrangeInterpolationNoBar      
!
!     ////////////////////////////////////////////////////////////////
!
!     -----------------------------------------------------------------------------
!!    Compute the value of the interpolant derivative using the barycentric formula
!     -----------------------------------------------------------------------------
!
      REAL(KIND=RP) FUNCTION LagrangeInterpolantDerivative( x, N, nodes, values, weights)
!
!     ---------
!     Arguments
!     ---------
!
      use Utilities, only: almostEqual
      REAL(KIND=RP)                 :: x
      INTEGER                       :: N
      REAL(KIND=RP), DIMENSION(0:N) :: nodes
      REAL(KIND=RP), DIMENSION(0:N) :: values
      REAL(KIND=RP), DIMENSION(0:N) :: weights
!
!     ---------------
!     Local Variables
!     ---------------
!
      INTEGER       :: j, i
      REAL(KIND=RP) :: t, numerator, denominator
      REAL(KIND=RP) :: p
      LOGICAL       :: atNode
!
!     --------------------------
!     See if the point is a node
!     --------------------------
!
      atNode = .FALSE.
      numerator   = 0.0_RP
      DO j = 0, N 
         IF( AlmostEqual( x, nodes(j) ) )    THEN
            atNode      = .TRUE.
            p           = values(j)
            denominator = -weights(j)
            i           = j
            EXIT 
         END IF
      END DO
      
      IF ( atNode )     THEN
         DO j = 0, N
            IF( j == i )    CYCLE
            numerator = numerator + weights(j)*(p - values(j))/( x - nodes(j) )
         END DO
      ELSE
         p = LagrangeInterpolation( x, N, nodes, values, weights)
         denominator = 0.0_RP
         DO j = 0, N
            t = weights(j)/( x - nodes(j) )
            numerator = numerator + t*(p - values(j))/( x - nodes(j) )
            denominator = denominator + t
         END DO
      END IF
      LagrangeInterpolantDerivative = numerator/denominator

      END FUNCTION LagrangeInterpolantDerivative
!
!     ////////////////////////////////////////////////////////////////
!
!     ------------------------------------------------------------------
!!    Compute the matrix T that transforms between the old set of nodes
!!    to the new set.
!     ------------------------------------------------------------------
!
      pure SUBROUTINE PolynomialInterpolationMatrix( N, M, oldNodes, weights, newNodes, T)
!
!     ---------
!     Arguments
!     ---------
!
      use Utilities, only: almostEqual
      INTEGER                      , INTENT(IN)  :: N, M
      REAL(KIND=RP), DIMENSION(0:N), INTENT(IN)  :: oldNodes
      REAL(KIND=RP), DIMENSION(0:M), INTENT(IN)  :: newNodes
      REAL(KIND=RP), DIMENSION(0:N), INTENT(IN)  :: weights
      
      REAL(KIND=RP), DIMENSION(0:M,0:N), INTENT(OUT) :: T
!
!     ---------------
!     Local Variables
!     ---------------
!
      INTEGER           :: j,k
      REAL(KIND=RP)     :: s, tmp
      LOGICAL           :: rowHasMatch
     
      DO k = 0,M
      
         rowHasMatch = .FALSE.
         DO j = 0,N
            T(k,j) = 0.0_RP
            IF( AlmostEqual(newNodes(k),oldNodes(j)) ) THEN
               rowHasMatch = .TRUE.
               T(k,j) = 1.0_RP
            END IF
         END DO 
         
         IF( .NOT.rowHasMatch )     THEN
            s = 0.0_RP
            DO j = 0,N
               tmp    = weights(j)/( newNodes(k) - oldNodes(j) )
               T(k,j) = tmp
               s      = s + tmp
            END DO
            DO j = 0, N 
               T(k,j) = T(k,j)/s
            END DO
         END IF
      END DO
      END SUBROUTINE PolynomialInterpolationMatrix
!
!////////////////////////////////////////////////////////////////////////
!
   SUBROUTINE Create1DInterpolationMatrix(Mat,N1,N2,x1,x2)
      IMPLICIT NONE
!
!     -----------------------------------------------------------
!     Creates a 3D Lagrange interpolation matrix from a grid with 
!     coordinates x1, y1, z1 (origin) to a grid with coordinates
!     x2, y2, z2 (destination)
!     -----------------------------------------------------------
!
      REAL(KIND=RP), ALLOCATABLE  :: Mat(:,:)     !>  Interpolation matrix
      INTEGER                     :: N1  !<  Origin order
      INTEGER                     :: N2  !<  Destination order
      REAL(KIND=RP), DIMENSION(:) :: x1(0:N1)     !<  Nodes in origin
      REAL(KIND=RP), DIMENSION(:) :: x2(0:N2)     !<  Nodes in destination
      !----------------------------------------------------------
      INTEGER :: i,j              ! Coordinate counters
      !----------------------------------------------------------
      
      ALLOCATE(Mat(N2 + 1,N1 + 1))
      
      DO j=0, N1                                  ! Column index   
         DO i=0, N2                               ! Row index
            Mat(i+1,j+1) =  LagrangeInterpolationNoBar(x2(i),N1,x1,j)
         END DO
      END DO
      
   END SUBROUTINE Create1DInterpolationMatrix
!
!////////////////////////////////////////////////////////////////////////
!
   SUBROUTINE Create1DRestrictionMatrix(Mat,N1,N2,x1,x2,w1,w2)
      IMPLICIT NONE
!
!     -----------------------------------------------------------
!     Creates a 3D Lagrange interpolation matrix from a grid with 
!     coordinates x1, y1, z1 (origin) to a grid with coordinates
!     x2, y2, z2 (destination)
!     -----------------------------------------------------------
!
      REAL(KIND=RP), ALLOCATABLE  :: Mat(:,:)     !>  Interpolation matrix
      INTEGER                     :: N1  !<  Origin order
      INTEGER                     :: N2  !<  Destination order
      REAL(KIND=RP), DIMENSION(:) :: x1(0:N1)     !<  Nodes in origin
      REAL(KIND=RP), DIMENSION(:) :: x2(0:N2)     !<  Nodes in destination
      REAL(KIND=RP), DIMENSION(:) :: w1(0:N1)     !<  Nodes in origin
      REAL(KIND=RP), DIMENSION(:) :: w2(0:N2)     !<  Nodes in destination
      !----------------------------------------------------------
      INTEGER :: i,j                    ! Coordinate counters
      !----------------------------------------------------------
      
      ALLOCATE(Mat(N2 + 1,N1 + 1))
      
      DO j=0, N1                                  ! Column index   
         DO i=0, N2                               ! Row index
            Mat(i+1,j+1) =  LagrangeInterpolationNoBar(x1(j),N2,x2,i) * w1(j)
         END DO
      END DO
      
      ! Create Mass matrix and finish computing interpolation operator
      DO i=0, N2                                  ! Row index
         ! Matrix Multiplication I = M⁻¹S (taking advantage of the diagonal matrix)
         Mat(i+1,:) = Mat(i+1,:) / w2(i)
      END DO
      
   END SUBROUTINE Create1DRestrictionMatrix
!
!////////////////////////////////////////////////////////////////////////
!
   SUBROUTINE Create3DInterpolationMatrix(Mat,N1x,N1y,N1z,N2x,N2y,N2z,x1,y1,z1,x2,y2,z2)
      IMPLICIT NONE
!
!     -----------------------------------------------------------
!     Creates a 3D Lagrange interpolation matrix from a grid with 
!     coordinates x1, y1, z1 (origin) to a grid with coordinates
!     x2, y2, z2 (destination)
!     -----------------------------------------------------------
!
      REAL(KIND=RP), ALLOCATABLE  :: Mat(:,:)     !>  Interpolation matrix
      INTEGER                     :: N1x,N1y,N1z  !<  Origin order
      INTEGER                     :: N2x,N2y,N2z  !<  Destination order
      REAL(KIND=RP), DIMENSION(:) :: x1(0:N1x),y1(0:N1y),z1(0:N1z)     !<  Nodes in origin
      REAL(KIND=RP), DIMENSION(:) :: x2(0:N2x),y2(0:N2y),z2(0:N2z)     !<  Nodes in destination
      !----------------------------------------------------------
      INTEGER :: i,j,k,l,m,n      ! Coordinate counters
      INTEGER :: r,s              ! Matrix index counters
      INTEGER :: NDOFEL1, NDOFEL2 ! Degrees of freedom in origin and destination
      !----------------------------------------------------------
      
      NDOFEL2 = (N2x + 1) * (N2y + 1) * (N2z + 1)
      NDOFEL1 = (N1x + 1) * (N1y + 1) * (N1z + 1)
      ALLOCATE(Mat(NDOFEL2,NDOFEL1))
      
      DO k=0, N1z
         DO j=0, N1y
            DO i=0, N1x
               r = i + j*(N1x + 1) + k*(N1x + 1)*(N1y + 1) + 1            ! Column index
               DO n=0, N2z
                  DO m=0, N2y
                     DO l=0, N2x
                        s = l + m*(N2x + 1) + n*(N2x + 1)*(N2y + 1) + 1   ! Row index
                        
                        Mat(s,r) =  LagrangeInterpolationNoBar(x2(l),N1x,x1,i) * &
                                    LagrangeInterpolationNoBar(y2(m),N1y,y1,j) * &
                                    LagrangeInterpolationNoBar(z2(n),N1z,z1,k)
                     END DO
                  END DO
               END DO
            END DO
         END DO
      END DO
   END SUBROUTINE Create3DInterpolationMatrix
!
!////////////////////////////////////////////////////////////////////////
!
   SUBROUTINE Create3DRestrictionMatrix(Mat,N1x,N1y,N1z,N2x,N2y,N2z,x1,y1,z1,x2,y2,z2,w1x,w1y,w1z,w2x,w2y,w2z)
      IMPLICIT NONE
!
!     -----------------------------------------------------------
!     Creates an L2-3D Lagrange interpolation matrix from a grid  
!     with coordinates x1, y1, z1 (origin) to a grid with coordinates
!     x2, y2, z2 (destination)
!     -----------------------------------------------------------
!
      REAL(KIND=RP), ALLOCATABLE  :: Mat(:,:)     !>  Interpolation matrix
      INTEGER                     :: N1x,N1y,N1z  !<  Origin order
      INTEGER                     :: N2x,N2y,N2z  !<  Destination order
      REAL(KIND=RP), DIMENSION(:) :: x1 (0:N1x),y1 (0:N1y),z1 (0:N1z)     !<  Nodes in origin
      REAL(KIND=RP), DIMENSION(:) :: x2 (0:N2x),y2 (0:N2y),z2 (0:N2z)     !<  Nodes in destination
      REAL(KIND=RP), DIMENSION(:) :: w1x(0:N1x),w1y(0:N1y),w1z(0:N1z)     !<  Weights in origin
      REAL(KIND=RP), DIMENSION(:) :: w2x(0:N2x),w2y(0:N2y),w2z(0:N2z)     !<  Weights in destination
      !----------------------------------------------------------
      INTEGER       :: i,j,k,l,m,n      ! Coordinate counters
      INTEGER       :: r,s              ! Matrix index counters
      INTEGER       :: NDOFEL1, NDOFEL2 ! Degrees of freedom in origin and destination
      REAL(KIND=RP) :: MASSterm         ! 
      !----------------------------------------------------------
      
      NDOFEL2 = (N2x + 1) * (N2y + 1) * (N2z + 1)
      NDOFEL1 = (N1x + 1) * (N1y + 1) * (N1z + 1)
      ALLOCATE(Mat(NDOFEL2,NDOFEL1))
      
      Mat = 0.0_RP
      
      ! Create S matrix and store it directly in "Mat"
      DO k=0, N1z
         DO j=0, N1y
            DO i=0, N1x
               r = i + j*(N1x + 1) + k*(N1x + 1)*(N1y + 1) + 1            ! Column index
               DO n=0, N2z
                  DO m=0, N2y
                     DO l=0, N2x
                        s = l + m*(N2x + 1) + n*(N2x + 1)*(N2y + 1) + 1   ! Row index
                        
                        Mat(s,r) = LagrangeInterpolationNoBar(x1(i),N2x,x2,l) * &
                                   LagrangeInterpolationNoBar(y1(j),N2y,y2,m) * &
                                   LagrangeInterpolationNoBar(z1(k),N2z,z2,n) * &
                                   w1x(i) * w1y(j) * w1z(k)
                     END DO
                  END DO
               END DO
            END DO
         END DO
      END DO
      
      ! Create Mass matrix and finish computing interpolation operator
      DO n=0, N2z
         DO m=0, N2y
            DO l=0, N2x
               s = l + m*(N2x + 1) + n*(N2x + 1)*(N2y + 1) + 1   ! Row index
      
               MASSterm = w2x(l) * w2y(m) * w2z(n)
               
               ! Matrix Multiplication I = M⁻¹S (taking advantage of the diagonal matrix)
               Mat(s,:) = Mat(s,:) / MASSterm
            END DO
         END DO
      END DO
      
   END SUBROUTINE Create3DRestrictionMatrix
!
!////////////////////////////////////////////////////////////////////////
!
   SUBROUTINE Interpolate3D(Q1, Q2, Interp, N1x, N1y, N1z, N2x, N2y, N2z)
      IMPLICIT NONE
      INTEGER        :: N1x, N1y, N1z
      INTEGER        :: N2x, N2y, N2z                      !<  Polynomial orders
      REAL(KIND=RP)  :: Q1((N1x+1)*(N1y+1)*(N1z+1))        !<  Solution to be interpolated (grid (1))
      REAL(KIND=RP)  :: Q2((N2x+1)*(N2y+1)*(N2z+1))        !>  Interpolated solution       (grid (2))
      REAL(KIND=RP)  :: Interp((N2x+1)*(N2y+1)*(N2z+1), & 
                               (N1x+1)*(N1y+1)*(N1z+1))    !<  Interpolation matrix
                               
      Q2 = MATMUL(Interp,Q1)
   END SUBROUTINE Interpolate3D
!
!////////////////////////////////////////////////////////////////////////
!
      SUBROUTINE InterpolateToNewPoints( N, M, T, f, fInterp )
!
!-------------------------------------------------------------------
!!    Use matrix Multiplication to interpolate between two sets of
!!    nodes.
!-------------------------------------------------------------------
!
!
!     ---------
!     Arguments
!     ---------
!
      INTEGER                          , INTENT(IN)  :: N, M
      REAL(KIND=RP), DIMENSION(0:M,0:N), INTENT(IN)  :: T
      REAL(KIND=RP), DIMENSION(0:N)    , INTENT(IN)  :: f
      REAL(KIND=RP), DIMENSION(0:M)    , INTENT(OUT) :: fInterp
!
!     ---------------
!     Local Variables
!     ---------------
!
      INTEGER       :: i,j
      REAL(KIND=RP) :: tmp
      
      DO i = 0, M
         tmp = 0.0_RP
         DO j = 0, N
            tmp = tmp + T(i,j)*f(j)
         END DO 
         fInterp(i) = tmp
      END DO 
      
      END SUBROUTINE InterpolateToNewPoints
!                                                                       
!///////////////////////////////////////////////////////////////////////
!
      SUBROUTINE Interp3DArray( inDim, inArray, outDim, outArray, interpXi, interpEta, interpZeta )
      

      IMPLICIT NONE
!
!     ---------
!     Arguments
!     ---------
!
      INTEGER          , DIMENSION(3)                              :: inDim
      INTEGER          , DIMENSION(3)                              :: outDim
      REAL(KIND=RP), DIMENSION(0:inDim(1) , 0:inDim(2) , 0:inDim(3) ) :: inArray
      REAL(KIND=RP), DIMENSION(0:outDim(1), 0:outDim(2), 0:outDim(3)) :: outArray
      REAL(KIND=RP), DIMENSION(0:outDim(1), 0:inDim(1))            :: interpXi
      REAL(KIND=RP), DIMENSION(0:outDim(2), 0:inDim(2))            :: interpEta
      REAL(KIND=RP), DIMENSION(0:outDim(3), 0:inDim(3))            :: interpZeta
!
!     ---------------
!     Local Variables
!     ---------------
!
      REAL(KIND=RP), DIMENSION(:)    , POINTER :: tempIn,tempOut
      REAL(KIND=RP), DIMENSION(:,:,:), POINTER :: tempArray
      
      INTEGER          , DIMENSION(3)              :: maxDim
      INTEGER                                      :: i,j,k
!
!     -------------------------------------------------
!     Allocate temporary storage for the interpolations
!     -------------------------------------------------
!
      ALLOCATE( tempIn ( 0:MAXVAL(inDim)  ) )
      ALLOCATE( tempOut( 0:MAXVAL(outDim) ) )
      maxDim = MAX(inDim,OutDim)
      ALLOCATE( tempArray(0:maxDim(1), 0:maxDim(2), 0:maxDim(3)) )
!
!     -----------------
!     Interpolate in Xi
!     -----------------
!
      DO k = 0, inDim(3)
         DO j = 0, inDim(2)
            DO i = 0, inDim(1)
               tempIn(i) = inArray(i,j,k)
            END DO 
            CALL InterpolateToNewPoints( inDim(1), outDim(1), interpXi, tempIn, tempOut )
            DO i = 0, outDim(1)
               tempArray(i,j,k) = tempOut(i)
            END DO 
         END DO
      END DO 
!
!     ------------------
!     Interpolate in Eta
!     ------------------
!
      DO k = 0, inDim(3)
         DO i = 0, outDim(1)
            DO j = 0, inDim(2)
               tempIn(j) = tempArray(i,j,k)
            END DO 
            CALL InterpolateToNewPoints( inDim(2), outDim(2), interpEta, tempIn, tempOut )
            DO j = 0, outDim(2)
               tempArray(i,j,k) = tempOut(j)
            END DO 
         END DO 
      END DO 
!
!     -------------------
!     Interpolate in Zeta
!     -------------------
!
      DO j = 0, outDim(2)
         DO i = 0, outDim(1)
            DO k = 0, inDim(3)
               tempIn(k) = tempArray(i,j,k)
            END DO 
            CALL InterpolateToNewPoints( inDim(3), outDim(3), interpZeta, tempIn, tempOut )
            DO k = 0, outDim(3)
               outArray(i,j,k) = tempOut(k)
            END DO 
         END DO 
      END DO 
!
!     ---------------
!     Clean up memory
!     ---------------
!
      DEALLOCATE (tempIn)
      DEALLOCATE (tempOut)
      DEALLOCATE (tempArray)
      
      END SUBROUTINE Interp3DArray      
!                                                                       
!///////////////////////////////////////////////////////////////////////
!

!
!    /////////////////////////////////////////////////////////////////
!
!     ----------------------------------------------------------------
!!    Compute the derivative matrix.
!!    D(i,j) = \prime \ell_j(x_i)
!     ----------------------------------------------------------------
!
      SUBROUTINE PolynomialDerivativeMatrix( N, nodes, D )
!
!     ---------
!     Arguments
!     ---------
!
      INTEGER                           , INTENT(IN)  :: N
      REAL(KIND=RP), DIMENSION(0:N)     , INTENT(IN)  :: nodes
      REAL(KIND=RP), DIMENSION(0:N, 0:N), INTENT(OUT) :: D
!
!     ---------------
!     Local Variables
!     ---------------
!
      INTEGER                        :: i, j
      REAL(KIND=RP), DIMENSION(0:N)  :: baryWeights
!      
      CALL BarycentricWeights( N, nodes, baryWeights )
      
      DO i = 0, N
         D(i,i) = 0.0_RP
         DO j = 0, N
            IF( j /= i )     THEN
               D(i,j) = baryWeights(j)/( baryWeights(i)*(nodes(i) - nodes(j)) )
               D(i,i) = D(i,i) - D(i,j)
            END IF
         END DO
      END DO

      END SUBROUTINE PolynomialDerivativeMatrix
!
!    /////////////////////////////////////////////////////////////////
!
!     ----------------------------------------------------------------
!!    Compute the (transpose) of the m^{th} derivative matrix.
!     ----------------------------------------------------------------
!
      SUBROUTINE mthPolynomialDerivativeMatrix( m, N, nodes, D )
!
!     ---------
!     Arguments
!     ---------
!
      INTEGER                           , INTENT(IN)  :: m
      INTEGER                           , INTENT(IN)  :: N
      REAL(KIND=RP), DIMENSION(0:N)     , INTENT(IN)  :: nodes
      REAL(KIND=RP), DIMENSION(0:N, 0:N), INTENT(OUT) :: D
!
!     ---------------
!     Local Variables
!     ---------------
!
      INTEGER                            :: i, j, l
      REAL(KIND=RP), DIMENSION(0:N)      :: baryWeights
      REAL(KIND=RP), DIMENSION(0:N, 0:N) :: DOld
!      
      CALL BarycentricWeights( N, nodes, baryWeights )
      CALL PolynomialDerivativeMatrix( N, nodes, DOld )
      
      DO l = 2, m
         DO i = 0, N
            D(i,i) = 0.0_RP
            DO j = 0, N
               IF( j /= i )     THEN
                  D(j,i) = l*(baryWeights(j)*DOld(i,i)/baryWeights(i) - DOld(j,i))/(nodes(i) - nodes(j))
                  D(i,i) = D(i,i) - D(j,i)
               END IF
            END DO
         END DO
         DOld = D
      END DO

      END SUBROUTINE mthPolynomialDerivativeMatrix
!
! /////////////////////////////////////////////////////////////////////
!
!---------------------------------------------------------------------
!!    Compute a square matrix times a vector where the vector is
!!    stored as the first dimension of a two-dimensional array.
!---------------------------------------------------------------------
!
      SUBROUTINE MMMultiply2D1(A,n,x,y)
      IMPLICIT NONE
!
!     ---------
!     Arguments
!     ---------
!
      INTEGER      , DIMENSION(2)            , INTENT(IN)  :: n
      REAL(KIND=RP), DIMENSION(0:n(1),0:n(1)), INTENT(IN)  :: A
      REAL(KIND=RP), DIMENSION(0:n(1),0:n(2)), INTENT(IN)  :: x
      
      REAL(KIND=RP), DIMENSION(0:n(1),0:n(2)), INTENT(OUT) :: y
!
!     ---------------
!     Local Variables
!     ---------------
!
      INTEGER           :: i,j,s
!
      y = 0.0_RP
      
      DO j = 0, n(2)
         DO s = 0, n(1)
            DO i = 0, n(1)
               y(i,j) = y(i,j) + A(i,s)*x(s,j)
            END DO
         END DO
      END DO
      
      END SUBROUTINE MMMultiply2D1
!
! /////////////////////////////////////////////////////////////////////
!
!---------------------------------------------------------------------
!!    Compute a square matrix times a vector where the vector is
!!    stored as the second dimension of a two-dimensional array.
!---------------------------------------------------------------------
!
      SUBROUTINE MMMultiply2D2(A,n,x,y)
      IMPLICIT NONE
!
!     ---------
!     Arguments
!     ---------
!
      INTEGER      , DIMENSION(2)            , INTENT(IN)  :: n
      REAL(KIND=RP), DIMENSION(0:n(2),0:n(2)), INTENT(IN)  :: A
      REAL(KIND=RP), DIMENSION(0:n(1),0:n(2)), INTENT(IN)  :: x
      
      REAL(KIND=RP), DIMENSION(0:n(1),0:n(2)), INTENT(OUT) :: y
!
!     ---------------
!     Local Variables
!     ---------------
!
      INTEGER           :: i,j,s
!
      y = 0.0_RP
      
      DO s = 0, n(2)
         DO j = 0, n(2)
            DO i = 0, n(1)
               y(i,j) = y(i,j) + A(j,s)*x(i,s)
            END DO
         END DO
      END DO
      
      END SUBROUTINE MMMultiply2D2
!
!////////////////////////////////////////////////////////////////////////////////////////
!
      SUBROUTINE PolyDirectMatrixMultiplyDeriv( f, fDeriv, DT, N )
!
!     Compute the derivative approximation by simple matrix multiplication
!     This routine assumes that the matrix has been transposed for faster
!     inner loop access.
!
!     ----------
!     Arguments:
!     ----------
!
      INTEGER                          , INTENT(IN)  :: N
      REAL(KIND=RP), DIMENSION(0:N)    , INTENT(IN)  :: f
      REAL(KIND=RP), DIMENSION(0:N,0:N), INTENT(IN)  :: DT

      REAL(KIND=RP), DIMENSION(0:N)    , INTENT(OUT) :: fDeriv
!
!     ----------------
!     Local Variables:
!     ----------------
!
      INTEGER       :: i, j
      REAL(KIND=RP) :: t
      
      DO i = 0, N
         t = 0.0_RP
         DO j = 0, N
            t = t + DT(j,i)*f(j)
         END DO
         fDeriv(i) = t
      END DO
      
      END SUBROUTINE PolyDirectMatrixMultiplyDeriv
!
!////////////////////////////////////////////////////////////////////////////////////////
!
      SUBROUTINE MatrixMultiplyDeriv( f, fDeriv, D, N, transp )
!
!     Compute the derivative approximation by simple matrix multiplication
!     This routine assumes that the matrix has been transposed for faster
!     inner loop access.
!
!     ----------
!     Arguments:
!     ----------
!
      INTEGER                          , INTENT(IN)  :: N
      REAL(KIND=RP), DIMENSION(0:N)    , INTENT(IN)  :: f
      REAL(KIND=RP), DIMENSION(0:N,0:N), INTENT(IN)  :: D
      INTEGER                          , INTENT(IN)  :: transp

      REAL(KIND=RP), DIMENSION(0:N)    , INTENT(OUT) :: fDeriv
!
!     ----------------
!     Local Variables:
!     ----------------
!
      INTEGER       :: i, j
      REAL(KIND=RP) :: t
      
      SELECT CASE ( transp )
      
         CASE( MXV_DIRECT )
             DO i = 0, N
               t = 0.0_RP
               DO j = 0, N
                  t = t + D(j,i)*f(j)
               END DO
               fDeriv(i) = t
            END DO
        
         CASE( MXV_TRANSPOSE )
             DO i = 0, N
               t = 0.0_RP
               DO j = 0, N
                  t = t + D(i,j)*f(j)
               END DO
               fDeriv(i) = t
            END DO
         
         CASE DEFAULT
      END SELECT
      
      END SUBROUTINE MatrixMultiplyDeriv
!
!////////////////////////////////////////////////////////////////////////////////////////
!
      SUBROUTINE EOMatrixDerivative( u, uDeriv, D, N )
      IMPLICIT NONE 
!
!     Compute the matrix derivative using the even-odd decomposition    
!     e = even terms                                                    
!     ep = even derivative                                              
!     o = odd terms                                                     
!     op = odd derivative
!
      INTEGER                          , INTENT(IN)  :: N
      REAL(KIND=RP), DIMENSION(0:N)    , INTENT(IN)  :: u
      REAL(KIND=RP), DIMENSION(0:N,0:N), INTENT(IN)  :: D
      
      REAL(KIND=RP), DIMENSION(0:N)    , INTENT(OUT) :: uDeriv

      REAL(KIND=RP), DIMENSION(0:N) :: e, o, ep, op
      REAL(KIND=RP)                 :: sume, sumo
      INTEGER                       :: M, nHalf, I, J
!
!     ----------------------------
!     Compute even and odd vectors                                      
!     ----------------------------
!
      M = (N+1)/2
      DO  j = 0,M
         e(j) = 0.5*(u(j) + u(n-j)) 
         o(j) = 0.5*(u(j) - u(n-j)) 
      END DO
!
!     -------------------------------------
!     Compute even and odd derivative terms                             
!     -------------------------------------
!
      DO i = 0, M-1
         sume = 0.0 
         sumo = 0.0 
         DO j = 0, M-1
            sume = sume + (D(j,i) + D(n-j,i))*e(j) 
            sumo = sumo + (D(j,i) - D(n-j,i))*o(j) 
         END DO
         ep(i) = sume 
         op(i) = sumo
      END DO
!
!     --------------------------------
!     Add in middle term if n+1 is odd                                   
!     --------------------------------
!
      nHalf = (n+1)/2
      IF(2*nHalf /= n+1)     THEN 
         DO i = 0, M-1 
            ep(i) = ep(i) + D(M,i)*e(M) 
         END DO
         ep(M) = 0.0 
         i = nHalf 
         sumo = 0.0 
         DO j = 0,M-1 
            sumo = sumo + (D(j,M) - D(n-j,M))*o(j) 
         END DO
         op(i) = sumo 
      END IF 
!
!     --------------------------
!     Construct full derivative                                         
!     --------------------------
!
      DO j = 0,M-1
         uDeriv(j)   = ( ep(j) + op(j)) 
         uDeriv(n-j) = (-ep(j) + op(j)) 
      END DO

      IF(2*nHalf /= n+1)     THEN 
         uDeriv(M) = (ep(M) + op(M))
      END IF 
!
      END SUBROUTINE EOMatrixDerivative                                           
!
!    /////////////////////////////////////////////////////////////////
!
!     ----------------------------------------------------------------
!     Discrete Galerkin(Continuous) Derivative matrix
!     for the second derivative.
!     ----------------------------------------------------------------
!
      SUBROUTINE DiscreteGalerkinDerivMatrix( N, G )
      USE GaussQuadrature
!
!     ---------
!     Arguments
!     ---------
!
      INTEGER                          , INTENT(IN)  :: N
      REAL(KIND=RP), DIMENSION(0:N,0:N), INTENT(OUT) :: G
!
!     ---------------
!     Local Variables
!     ---------------
!
     INTEGER                           :: k,m,j
     REAL(KIND=RP), DIMENSION(0:N)     :: x, w
     REAL(KIND=RP), DIMENSION(0:N,0:N) :: D
     REAL(KIND=RP)                     :: s
!      
     CALL LegendreLobattoNodesAndWeights( N, x, w )
     CALL PolynomialDerivativeMatrix( N, x, D )
     
     DO j = 0, N
        DO m = 0, N
           s = 0.0_RP
           DO k = 0, N
              s = s + D(m,k)*D(j,k)*w(k)
           END DO 
           G(m,j) = s/w(j)
        END DO 
     END DO 

      END SUBROUTINE DiscreteGalerkinDerivMatrix
!
! /////////////////////////////////////////////////////////////////////
!
   function JacobiPolynomial( N, x, nx, alpha, beta ) result (PolMat) 
!  ---------------------------------------------------------
!  General routine to evaluate orthonormal Jacobi polynomial
!  of type (alpha,beta) recursively at points defined by an
!  array x - essentially computes 1D Vandermonde matrix.   
!  ---------------------------------------------------------
      implicit none
!-----Input---------------------------------------------------------------
      integer,                            intent(in) :: N      !<  Polynomial order
      real(kind=RP), dimension(0:nx),     intent(in) :: x      !<  Where to evaluate the derivative
      integer,                            intent(in) :: nx 
      real(kind=RP),                      intent(in) :: alpha 
      real(kind=RP),                      intent(in) :: beta 
!-----Output--------------------------------------------------------------
      real(kind=RP), dimension(0:nx,0:N) :: PolMat !<  Resulting matrix with polynomials
!-----Local-Variables-----------------------------------------------------
      integer                        :: i,j
      real(kind=RP)                  :: gamma0, gamma1 
      real(kind=RP)                  :: aold, anew, bnew, tmp 
      real(kind=RP), dimension(0:nx) :: tmp2
!  ---------------------------------------------------------

      PolMat = 0.0_RP

      ! Compute initial value P1
      gamma0 = 2.0_RP**(alpha+beta+1)/(alpha+beta+1)*gamma(alpha+1)*gamma(beta+1)/gamma(alpha+beta+1)
      PolMat(:,0) = 1.0_RP/sqrt(gamma0)

      ! Compute initial value P2
      gamma1 = (alpha+1)*(beta+1)/(alpha+beta+3)*gamma0
      PolMat(:,1) = ((alpha+beta+2)*x/2 + (alpha-beta)/2)/sqrt(gamma1)

      ! Recurrent formula for the rest of polynomials
      aold = 2/(2+alpha+beta)*sqrt((alpha+1)*(beta+1)/(alpha+beta+3))
      do i=2,N 
         tmp = 2*(i-1)+alpha+beta
         anew = 2/(tmp+2)*sqrt( (i)*(i+alpha+beta)*(i+alpha)*(i+beta)/(tmp+1)/(tmp+3))
         bnew = - (alpha**2-beta**2)/tmp/(tmp+2)
         do j=0,nx
            tmp2(j) = (x(j)-bnew) * PolMat(j,i-1)
         end do
         PolMat(:,i) = 1/anew*( -aold*PolMat(:,i-2) + tmp2)
         aold = anew
      end do

   end function JacobiPolynomial
!
! /////////////////////////////////////////////////////////////////////
!
   function ComputeVandermonde( N ) result (V) 
!  ---------------------------------------------------------
!  Routine to compute a square Vandermonde matrix for 3D hexa.
!  ---------------------------------------------------------
      use GaussQuadrature
      implicit none
!-----Input---------------------------------------------------------------
      integer,       dimension(3),        intent(in)    :: N      !  Polynomial order
!-----Output--------------------------------------------------------------
      real(kind=RP), dimension(0:(N(1)+1)*(N(2)+1)*(N(3)+1)-1,0:(N(1)+1)*(N(2)+1)*(N(3)+1)-1) :: V ! Vandermonde matrix 
!-----Local-Variables-----------------------------------------------------
      integer :: i,j,k,l,Vsize 
      real(kind=RP), dimension(0:(N(1)+1)*(N(2)+1)*(N(3)+1)-1,0:N(1)) :: Vx ! 1D Vandermonde matrix in x 
      real(kind=RP), dimension(0:(N(1)+1)*(N(2)+1)*(N(3)+1)-1,0:N(2)) :: Vy ! 1D Vandermonde matrix in y 
      real(kind=RP), dimension(0:(N(1)+1)*(N(2)+1)*(N(3)+1)-1,0:N(3)) :: Vz ! 1D Vandermonde matrix in z 
      !-1D-GLL-quadrature-points-and-weights:---------
      real(kind=RP), dimension(0:N(1))                                :: x1d, wx 
      real(kind=RP), dimension(0:N(2))                                :: y1d, wy 
      real(kind=RP), dimension(0:N(3))                                :: z1d, wz 
      real(kind=RP), dimension(0:N(1),0:N(2),0:N(3))                  :: x3d,y3d,z3d ! 3d arrays for GLL points
      real(kind=RP), dimension(0:(N(1)+1)*(N(2)+1)*(N(3)+1)-1)        :: x,y,z ! 3d arrays reshaped 
!  ---------------------------------------------------------

      ! Compute 1D quadrature points and weights
      call LegendreLobattoNodesAndWeights( N(1), x1d, wx )
      call LegendreLobattoNodesAndWeights( N(2), y1d, wy )
      call LegendreLobattoNodesAndWeights( N(3), z1d, wz )

      ! Create 3D grids
      do i=0,N(3)
         x3d(:,:,i) = spread( x1d, 1, N(2) )
         y3d(:,:,i) = spread( y1d, 2, N(1) )
      end do
      do i=0,N(1)
         z3d(i,:,:) = spread( z1d, 1, N(2))
      end do 

      ! Reshape arrays
      x = reshape(x3d,(/ (N(1)+1)*(N(2)+1)*(N(3)+1) /)) 
      y = reshape(y3d,(/ (N(1)+1)*(N(2)+1)*(N(3)+1) /)) 
      z = reshape(z3d,(/ (N(1)+1)*(N(2)+1)*(N(3)+1) /)) 

      ! Compute 1D polynomials
      Vx = JacobiPolynomial(N(1),x,(N(1)+1)*(N(2)+1)*(N(3)+1)-1,0.0_RP,0.0_RP)
      Vy = JacobiPolynomial(N(2),y,(N(1)+1)*(N(2)+1)*(N(3)+1)-1,0.0_RP,0.0_RP)
      Vz = JacobiPolynomial(N(3),z,(N(1)+1)*(N(2)+1)*(N(3)+1)-1,0.0_RP,0.0_RP)

      ! Build 3D hexahedra Vandermonde matrix
      V = 0.0_RP 
      Vsize = 0;
      do i=0,N(1); do j=0,N(2); do k=0,N(3) 
         do l=0,size(V,1)-1
            V(l,Vsize) = Vx(l,i) * Vy(l,j) * Vz(l,k)
         end do
         Vsize = Vsize+1;
      end do; end do; end do

   end function ComputeVandermonde
!
!////////////////////////////////////////////////////////////////////////////////////////
!
   subroutine ComputeModalForm( u, uhat, Neqn, Nel, Nx, Ny, Nz )
!  ---------------------------------------------------------
!  Routine to compute modal form of any 1D array.
!  ---------------------------------------------------------
      use DenseMatUtilities
      implicit none
!-----Arguments-----------------------------------------------------------
      real(kind=RP), dimension(Nel*Neqn*(Nx+1)*(Ny+1)*(Nz+1)) , intent(in)    :: u    ! 1D array of nodal storage
      real(kind=RP), dimension(Nel*Neqn*(Nx+1)*(Ny+1)*(Nz+1)) , intent(inout) :: uhat ! 1D array of modal storage
      integer                                                 , intent(in)    :: Neqn ! no. equations
      integer                                                 , intent(in)    :: Nel  ! no. elements
      !-polynomial-orders-in-x/y/z:-------------------
      integer                                                 , intent(in)    :: Nx      
      integer                                                 , intent(in)    :: Ny
      integer                                                 , intent(in)    :: Nz
!-----Local-Variables-----------------------------------------------------
      real(kind=rp), dimension(0:(Nx+1)*(Ny+1)*(Nz+1)-1,0:(Nx+1)*(Ny+1)*(Nz+1)-1) :: V 
      integer                                                                     :: i,j,k,l,iel
      integer                                                                     :: elsize
!  ---------------------------------------------------------

      ! compute Vandermonde matrix on 3D hexahedral
      V = ComputeVandermonde( (/ Nx, Ny, Nz /) )  

      ! save Vandermonde matrix in the matrix form
      open(25, file = 'Vandermonde3DHex.txt')  
      do i=0,size(V,1)-1
            do j=0,size(V,2)-1
                  ! write(25,'(1E19.11)', advance="no") V(i,j) ! better precision 
                  write(25,'(1F8.3)', advance="no") V(i,j) 
            end do
            write(25,*)
      end do
      close(25) 

      ! compute V^(-1)
      V = inverse(V)

      elsize = (Nx+1)*(Ny+1)*(Nz+1)
      uhat = 0.0_RP
      do iEl = 1, Nel; do i=1,Neqn*elsize,Neqn; do j=0,Neqn-1
         uhat( 1 + j + (iEl-1) * (elsize*Neqn) : iEl * (elsize*Neqn) : Neqn ) = &
            matmul( V , u( 1 + j + (iEl-1) * (elsize*Neqn) : iEl * (elsize*Neqn) : Neqn ) )
      end do; end do; end do


      ! save uhat in the form of (Nel*elsize,Neqn) matrix
      open(26, file = 'Qhat.txt')  
      do i=1,Neqn*Nel*elsize,Neqn; do j=0,Neqn-1
            write(26,'(1E22.14)', advance="no") uhat(i+j)
      end do; write(26,*); end do
      close(26) 

   end subroutine ComputeModalForm
!
! /////////////////////////////////////////////////////////////////////
!
!
!
!  **********     
   END MODULE PolynomialInterpAndDerivsModule
!  **********     