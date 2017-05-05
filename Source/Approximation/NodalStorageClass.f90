!
!////////////////////////////////////////////////////////////////////////
!
!      NodalStorage.f95
!      Created: 2008-01-15 10:35:59 -0500 
!      By: David Kopriva 
!
!     Algorithms:
!        Algorithm 
!
!////////////////////////////////////////////////////////////////////////
!
      MODULE NodalStorageClass
      USE SMConstants
      USE PolynomialInterpAndDerivsModule
      USE GaussQuadrature
      IMPLICIT NONE 
!
!     -----
!     Class
!     -----
!
      TYPE NodalStorage
         INTEGER                                    :: N
         REAL(KIND=RP)                , ALLOCATABLE :: xi(:), eta(:), zeta(:), w(:)
         REAL(KIND=RP)                , ALLOCATABLE :: standardDerivativeMatrix(:,:), D(:,:)
         REAL(KIND=RP), DIMENSION(:,:), ALLOCATABLE :: v                      ! Interpolation vector
         REAL(KIND=RP), DIMENSION(:,:), ALLOCATABLE :: b                      ! Boundary vector
         LOGICAL                                    :: Constructed = .FALSE.  ! Has this combination already been constructed?
         
         CONTAINS
         PROCEDURE :: construct => ConstructNodalStorage
         PROCEDURE :: destruct  => DestructNodalStorage
         
      END TYPE NodalStorage
!      
!     ========
      CONTAINS 
!     ========
!
!////////////////////////////////////////////////////////////////////////
!
      SUBROUTINE ConstructNodalStorage( this, N )
      IMPLICIT NONE
      CLASS(NodalStorage)      :: this
      INTEGER                  :: N
      
      INTEGER       :: i,j
      REAL(KIND=RP) :: D(0:N,0:N)
      REAL(KIND=RP) :: wb(0:N)
      INTEGER, PARAMETER :: LEFT = 1, RIGHT = 2, TOP = 2, BOTTOM = 1
            
      this % N = N
      ALLOCATE( this % xi(0:N), this % eta(0:N), this % zeta(0:N) )
      ALLOCATE( this % w(0:N)     )
      ALLOCATE( this % D(0:N,0:N) )
      ALLOCATE( this % standardDerivativeMatrix(0:N,0:N) )
      ALLOCATE( this % v(0:N,2) )
      ALLOCATE( this % b(0:N,2) )
!
!     -----------------
!     Nodes and weights
!     -----------------
!
      CALL GaussLegendreNodesAndWeights( N, this % xi , this % w )
      this % eta  = this % xi
      this % zeta = this % xi
!
!     -----------------
!     Derivative Matrix
!     -----------------
!
      CALL PolynomialDerivativeMatrix( N, this % xi, D )
      this % standardDerivativeMatrix = D
      DO j = 0, N 
         DO i = 0, N 
            this % D(j,i) = -D(i,j)*this % w(j)/this % w(i)
         END DO
      END DO
!
!     ---------------------
!     Interpolation vectors
!     ---------------------
!
      CALL BarycentricWeights( N, this % xi, wb )
      
      CALL InterpolatingPolynomialVector(  1.0_RP, N, this % xi, wb, this % b(:,RIGHT) )
      CALL InterpolatingPolynomialVector( -1.0_RP, N, this % xi, wb, this % b(:,LEFT) )
      
      this % v = this % b
      
      this % b(0:N,LEFT)  = this % b(0:N,LEFT)/this % w
      this % b(0:N,RIGHT) = this % b(0:N,RIGHT)/this % w
      
      this % Constructed = .TRUE.
      
      END SUBROUTINE ConstructNodalStorage
!
!////////////////////////////////////////////////////////////////////////
!
      SUBROUTINE DestructNodalStorage( this )
         IMPLICIT NONE
         CLASS(NodalStorage) :: this
         DEALLOCATE( this % xi, this % eta, this % zeta )
         DEALLOCATE( this % w )
         DEALLOCATE( this % D, this % v, this % b )
         deallocate( this % standardDerivativeMatrix)
         this % Constructed = .FALSE.
      END SUBROUTINE DestructNodalStorage
      
      END Module NodalStorageClass
