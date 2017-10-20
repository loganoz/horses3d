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
         INTEGER                                    :: Nx,Ny,Nz                        ! Polynomial orders in every direction
         REAL(KIND=RP)                , ALLOCATABLE :: xi(:), eta(:), zeta(:)          ! Node position in every direction
         REAL(KIND=RP)                , ALLOCATABLE :: wx(:), wy(:), wz(:)             ! Weights in every direction
         REAL(KIND=RP)                , ALLOCATABLE :: wbx(:), wby(:), wbz(:)          ! Barycentric weights in every direction
         REAL(KIND=RP)                , ALLOCATABLE :: standardDerivativeMatrix(:,:)   ! Standard derivative matrix for the x direction (to deprecate?... This is only used in one of the Components' test cases)
         REAL(KIND=RP), DIMENSION(:,:), ALLOCATABLE :: vx, vy, vz                      ! Interpolation vector
         REAL(KIND=RP)                , ALLOCATABLE :: Dx(:,:), Dy(:,:), Dz(:,:)       ! DG derivative matrices in every direction
         REAL(KIND=RP), DIMENSION(:,:), ALLOCATABLE :: bx, by, bz                      ! Boundary vector
         LOGICAL                                    :: Constructed = .FALSE.           ! Has this nodal storage already been constructed?
         REAL(KIND=RP)                , ALLOCATABLE :: hatDx(:,:) , hatDy(:,:) , hatDz(:,:)
         
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
      SUBROUTINE ConstructNodalStorage( this, Nx, Ny, Nz )
      IMPLICIT NONE
      CLASS(NodalStorage)      :: this       !<> Nodal storage being constructed
      INTEGER                  :: Nx, Ny, Nz !<  Polynomial orders in the different directions
      !--------------------------------------
      INTEGER       :: i,j
      REAL(KIND=RP) :: Dx(0:Nx,0:Nx),Dy(0:Ny,0:Ny),Dz(0:Nz,0:Nz)
      REAL(KIND=RP) :: wbx(0:Nx), wby(0:Ny), wbz(0:Nz)
      INTEGER, PARAMETER :: LEFT = 1, RIGHT = 2, TOP = 2, BOTTOM = 1
      !--------------------------------------
      
      this % Nx = Nx
      this % Ny = Ny
      this % Nz = Nz
      
      ALLOCATE( this % xi(0:Nx), this % eta(0:Ny), this % zeta(0:Nz) )
      ALLOCATE( this % wx(0:Nx), this % wy (0:Ny), this % wz  (0:Nz) )
      ALLOCATE( this % wbx(0:Nx), this % wby (0:Ny), this % wbz  (0:Nz) )
      ALLOCATE( this % Dx(0:Nx,0:Nx), this % Dy(0:Ny,0:Ny), this % Dz(0:Nz,0:Nz) )
      ALLOCATE( this % standardDerivativeMatrix(0:Nx,0:Nx) )
      ALLOCATE( this % vx(0:Nx,2), this % vy(0:Ny,2), this % vz(0:Nz,2) )
      ALLOCATE( this % bx(0:Nx,2), this % by(0:Ny,2), this % bz(0:Nz,2) )
      ALLOCATE( this % hatDx ( 0:Nx,0:Nx ) , this % hatDy(0:Ny,0:Ny) , this % hatDz(0:Nz,0:Nz)  ) 
!
!     -----------------
!     Nodes and weights
!     -----------------
!
      CALL GaussLegendreNodesAndWeights( Nx, this % xi  , this % wx )
      CALL GaussLegendreNodesAndWeights( Ny, this % eta , this % wy )
      CALL GaussLegendreNodesAndWeights( Nz, this % zeta, this % wz )
!
!     -----------------
!     Derivative Matrix
!     -----------------
!
      ! x Direction
      CALL PolynomialDerivativeMatrix( Nx, this % xi, Dx )
      this % standardDerivativeMatrix = Dx
      DO j = 0, Nx 
         DO i = 0, Nx 
            this % Dx(j,i) = -Dx(i,j)*this % wx(j)/this % wx(i)
         END DO
      END DO
      
      ! y Direction
      CALL PolynomialDerivativeMatrix( Ny, this % eta, Dy )
      DO j = 0, Ny
         DO i = 0, Ny 
            this % Dy(j,i) = -Dy(i,j)*this % wy(j)/this % wy(i)
         END DO
      END DO
      
      ! z direction
      CALL PolynomialDerivativeMatrix( Nz, this % zeta, Dz )
      DO j = 0, Nz
         DO i = 0, Nz 
            this % Dz(j,i) = -Dz(i,j)*this % wz(j)/this % wz(i)
         END DO
      END DO

      this % hatDx = - this % Dx
      this % hatDy = - this % Dy
      this % hatDz = - this % Dz
!
!     ---------------------
!     Interpolation vectors
!     ---------------------
!
      ! x Direction
      CALL BarycentricWeights( Nx, this % xi, wbx )
      
      CALL InterpolatingPolynomialVector(  1.0_RP, Nx, this % xi, wbx, this % bx(:,RIGHT) )
      CALL InterpolatingPolynomialVector( -1.0_RP, Nx, this % xi, wbx, this % bx(:,LEFT)  )
      
      this % wbx = wbx
      this % vx = this % bx
      
      this % bx(0:Nx,LEFT)  = this % bx(0:Nx,LEFT) /this % wx
      this % bx(0:Nx,RIGHT) = this % bx(0:Nx,RIGHT)/this % wx
      
      ! y Direction
      CALL BarycentricWeights( Ny, this % eta, wby )
      
      CALL InterpolatingPolynomialVector(  1.0_RP, Ny, this % eta, wby, this % by(:,RIGHT) )
      CALL InterpolatingPolynomialVector( -1.0_RP, Ny, this % eta, wby, this % by(:,LEFT)  )
      
      this % wby = wby
      this % vy = this % by
      
      this % by(0:Ny,LEFT)  = this % by(0:Ny,LEFT) /this % wy
      this % by(0:Ny,RIGHT) = this % by(0:Ny,RIGHT)/this % wy
      
      ! z Direction
      CALL BarycentricWeights( Nz, this % zeta, wbz )
      
      CALL InterpolatingPolynomialVector(  1.0_RP, Nz, this % zeta, wbz, this % bz(:,RIGHT) )
      CALL InterpolatingPolynomialVector( -1.0_RP, Nz, this % zeta, wbz, this % bz(:,LEFT)  )
      
      this % wbz = wbz
      this % vz = this % bz
      
      this % bz(0:Nz,LEFT)  = this % bz(0:Nz,LEFT) /this % wz
      this % bz(0:Nz,RIGHT) = this % bz(0:Nz,RIGHT)/this % wz
!
!     ---------
!     All done!
!     ---------
!
      
      this % Constructed = .TRUE.
      
      END SUBROUTINE ConstructNodalStorage
!
!////////////////////////////////////////////////////////////////////////
!
      SUBROUTINE DestructNodalStorage( this )
         IMPLICIT NONE
         CLASS(NodalStorage) :: this
         DEALLOCATE( this % xi, this % eta, this % zeta )
         DEALLOCATE( this % wx, this % wy , this % wz  )
         DEALLOCATE( this % Dx, this % Dy , this % Dz  )
         DEALLOCATE( this % vx, this % vy , this % vz  )        
         DEALLOCATE( this % bx, this % by , this % bz  )
         DEALLOCATE( this % hatDx , this % hatDy , this % hatDz )
         deallocate( this % standardDerivativeMatrix)
         this % Constructed = .FALSE.
      END SUBROUTINE DestructNodalStorage
      
      END Module NodalStorageClass
