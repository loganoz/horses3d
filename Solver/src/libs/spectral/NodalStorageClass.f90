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
#include "Includes.h"
MODULE NodalStorageClass
   USE SMConstants
   USE PolynomialInterpAndDerivsModule
   USE GaussQuadrature
   IMPLICIT NONE 

   private
   public GAUSS, GAUSSLOBATTO, NodalStorage

   integer, parameter      :: GAUSS = 1
   integer, parameter      :: GAUSSLOBATTO = 2
!
!  -------------------
!  Nodal storage class
!  -------------------
!
   type NodalStorage
      logical                                    :: Constructed = .FALSE.     ! Constructed flag
      integer                                    :: nodes                     ! Either GAUSS or GAUSSLOBATTO
      integer                                    :: Nx,Ny,Nz                  ! Polynomial orders in every direction
      real(kind=RP), dimension(:)  , allocatable :: xi, eta, zeta             ! Node position in every direction
      real(kind=RP), dimension(:)  , allocatable :: wx, wy, wz                ! Weights in every direction
      real(kind=RP), dimension(:)  , allocatable :: wbx, wby, wbz             ! Barycentric weights in every direction
      real(kind=RP), dimension(:,:), allocatable :: vx, vy, vz                ! Interpolation vector
      real(kind=RP), dimension(:,:), allocatable :: Dx, Dy, Dz                ! DG derivative matrices in every direction
      real(kind=RP), dimension(:,:), allocatable :: DTx, DTy, DTz             ! Trasposed DG derivative matrices in every direction
      real(kind=RP), dimension(:,:), allocatable :: hatDx, hatDy, hatDz       ! Weak form derivative matrices
      real(kind=RP), dimension(:,:), allocatable :: sharpDx, sharpDy, sharpDz ! (Two times) the strong form derivative matrix
      real(kind=RP), dimension(:,:), allocatable :: bx, by, bz                ! Boundary vector
      contains
         procedure :: construct => ConstructNodalStorage
         procedure :: destruct  => DestructNodalStorage
         procedure :: lxi       => NodalStorage_getlxi
         procedure :: leta      => NodalStorage_getleta
         procedure :: lzeta     => NodalStorage_getlzeta
         procedure :: dlxi      => NodalStorage_getdlxi
         procedure :: dleta     => NodalStorage_getdleta
         procedure :: dlzeta    => NodalStorage_getdlzeta

   END TYPE NodalStorage
!      
!     ========
      CONTAINS 
!     ========
!
!////////////////////////////////////////////////////////////////////////
!
      SUBROUTINE ConstructNodalStorage( this, nodes, Nx, Ny, Nz )
      IMPLICIT NONE
      CLASS(NodalStorage)      :: this       !<> Nodal storage being constructed
      integer, intent(in)      :: nodes
      INTEGER, intent(in)      :: Nx, Ny, Nz !<  Polynomial orders in the different directions
      !--------------------------------------
      INTEGER       :: i,j
      REAL(KIND=RP) :: wbx(0:Nx), wby(0:Ny), wbz(0:Nz)
      INTEGER, PARAMETER :: LEFT = 1, RIGHT = 2, TOP = 2, BOTTOM = 1
      !--------------------------------------
!
!     Set the nodal storage as constructed
!     ------------------------------------      
      this % Constructed = .TRUE.

      this % nodes = nodes
      this % Nx = Nx
      this % Ny = Ny
      this % Nz = Nz
      
      ALLOCATE( this % xi(0:Nx), this % eta(0:Ny), this % zeta(0:Nz) )
      ALLOCATE( this % wx(0:Nx), this % wy (0:Ny), this % wz  (0:Nz) )
      ALLOCATE( this % wbx(0:Nx), this % wby (0:Ny), this % wbz  (0:Nz) )
      ALLOCATE( this % Dx(0:Nx,0:Nx), this % Dy(0:Ny,0:Ny), this % Dz(0:Nz,0:Nz) )
      ALLOCATE( this % DTx(0:Nx,0:Nx), this % DTy(0:Ny,0:Ny), this % DTz(0:Nz,0:Nz) )
      ALLOCATE( this % vx(0:Nx,2), this % vy(0:Ny,2), this % vz(0:Nz,2) )
      ALLOCATE( this % bx(0:Nx,2), this % by(0:Ny,2), this % bz(0:Nz,2) )
      ALLOCATE( this % hatDx ( 0:Nx,0:Nx ) , this % hatDy(0:Ny,0:Ny) , this % hatDz(0:Nz,0:Nz)  ) 
!
!     -----------------
!     Nodes and weights
!     -----------------
!
      select case (this % nodes)
      case (GAUSS)
         CALL GaussLegendreNodesAndWeights( Nx, this % xi  , this % wx )
         CALL GaussLegendreNodesAndWeights( Ny, this % eta , this % wy )
         CALL GaussLegendreNodesAndWeights( Nz, this % zeta, this % wz )
      case (GAUSSLOBATTO)
         CALL LegendreLobattoNodesAndWeights( Nx, this % xi  , this % wx )
         CALL LegendreLobattoNodesAndWeights( Ny, this % eta , this % wy )
         CALL LegendreLobattoNodesAndWeights( Nz, this % zeta, this % wz )
         allocate( this % sharpDx(0:Nx,0:Nx), this % sharpDy(0:Ny,0:Ny), this % sharpDz(0:Nz,0:Nz) )
      case default
         print*, "Undefined nodes choice"
         errorMessage(STD_OUT)
         stop
      end select
!
!     -----------------
!     Derivative Matrix
!     -----------------
!
      ! x Direction
      CALL PolynomialDerivativeMatrix( Nx, this % xi, this % Dx )
      DO j = 0, Nx 
         DO i = 0, Nx 
            this % DTx(i,j) = this % Dx(j,i)
            this % hatDx(i,j) = this % DTx(i,j) * this % wx(j) / this % wx(i)
         END DO
      END DO

      ! y Direction
      CALL PolynomialDerivativeMatrix( Ny, this % eta, this % Dy )
      DO j = 0, Ny 
         DO i = 0, Ny 
            this % DTy(i,j) = this % Dy(j,i)
            this % hatDy(i,j) = this % DTy(i,j) * this % wy(j) / this % wy(i)
         END DO
      END DO

      ! z Direction
      CALL PolynomialDerivativeMatrix( Nz, this % zeta, this % Dz )
      DO j = 0, Nz 
         DO i = 0, Nz 
            this % DTz(i,j) = this % Dz(j,i)
            this % hatDz(i,j) = this % DTz(i,j) * this % wz(j) / this % wz(i)
         END DO
      END DO
!
!     --------------------------------------------------------------
!     Construct the strong form derivative matrices (Skew-Symmetric)
!     --------------------------------------------------------------
!
      if ( this % nodes .eq. GAUSSLOBATTO ) then
         this % sharpDx = 2.0_RP * this % Dx
         this % sharpDx(0,0) = 2.0_RP * this % Dx(0,0) + 1.0_RP / this % wx(0)
         this % sharpDx(Nx,Nx) = 2.0_RP * this % Dx(Nx,Nx) - 1.0_RP / this % wx(Nx)

         this % sharpDy = 2.0_RP * this % Dy
         this % sharpDy(0,0) = 2.0_RP * this % Dy(0,0) + 1.0_RP / this % wy(0)
         this % sharpDy(Ny,Ny) = 2.0_RP * this % Dy(Ny,Ny) - 1.0_RP / this % wy(Ny)

         this % sharpDz = 2.0_RP * this % Dz
         this % sharpDz(0,0) = 2.0_RP * this % Dz(0,0) + 1.0_RP / this % wz(0)
         this % sharpDz(Nz,Nz) = 2.0_RP * this % Dz(Nz,Nz) - 1.0_RP / this % wz(Nz)
      end if
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

      END SUBROUTINE ConstructNodalStorage
!
!////////////////////////////////////////////////////////////////////////
!
      SUBROUTINE DestructNodalStorage( this )
         IMPLICIT NONE
         CLASS(NodalStorage) :: this
!
!        Attempting to destruct a non-constructed nodal storage
!        ------------------------------------------------------
         if (.not. this % constructed ) return
!
!        Destruct otherwise
!        ------------------
         this % constructed = .FALSE.
         DEALLOCATE( this % xi, this % eta, this % zeta )
         DEALLOCATE( this % wx, this % wy , this % wz  )
         DEALLOCATE( this % Dx, this % Dy , this % Dz  )
         DEALLOCATE( this % DTx, this % DTy , this % DTz  )
         DEALLOCATE( this % hatDx , this % hatDy , this % hatDz )
         DEALLOCATE( this % vx, this % vy , this % vz  )        
         DEALLOCATE( this % bx, this % by , this % bz  )
         safedeallocate( this % sharpDx )    !  This matrices are just generated
         safedeallocate( this % sharpDy )    !  for Gauss-Lobatto discretizations.
         safedeallocate( this % sharpDz )    !

      END SUBROUTINE DestructNodalStorage
!
!////////////////////////////////////////////////////////////////////////
!
      function NodalStorage_getlxi(self, xi)
         implicit none
         class(NodalStorage), intent(in)  :: self
         real(kind=RP),       intent(in)  :: xi
         real(kind=RP)                    :: NodalStorage_getlxi(0:self % Nx)

         call InterpolatingPolynomialVector(xi , self % Nx , self % xi , self % wbx , NodalStorage_getlxi)

      end function NodalStorage_getlxi

      function NodalStorage_getleta(self, eta)
         implicit none
         class(NodalStorage), intent(in)  :: self
         real(kind=RP),       intent(in)  :: eta
         real(kind=RP)                    :: NodalStorage_getleta(0:self % Ny)

         call InterpolatingPolynomialVector(eta , self % Ny , self % eta , self % wby , NodalStorage_getleta)

      end function NodalStorage_getleta

      function NodalStorage_getlzeta(self, zeta)
         implicit none
         class(NodalStorage), intent(in)  :: self
         real(kind=RP),       intent(in)  :: zeta
         real(kind=RP)                    :: NodalStorage_getlzeta(0:self % Nz)

         call InterpolatingPolynomialVector(zeta , self % Nz , self % zeta , self % wbz , NodalStorage_getlzeta)

      end function NodalStorage_getlzeta

      function NodalStorage_getdlxi(self, xi)
         implicit none
         class(NodalStorage), intent(in)  :: self
         real(kind=RP),       intent(in)  :: xi
         real(kind=RP)                    :: NodalStorage_getdlxi(0:self % Nx)
!
!        ---------------
!        Local variables
!        ---------------
!
         integer  :: i

         do i = 0 , self % Nx
            NodalStorage_getdlxi(i) = EvaluateLagrangePolyDerivative( i, xi, self % Nx , self % xi)
         end do

      end function NodalStorage_getdlxi

      function NodalStorage_getdleta(self, eta)
         implicit none
         class(NodalStorage), intent(in)  :: self
         real(kind=RP),       intent(in)  :: eta
         real(kind=RP)                    :: NodalStorage_getdleta(0:self % Ny)
!
!        ---------------
!        Local variables
!        ---------------
!
         integer  :: i

         do i = 0 , self % Ny
            NodalStorage_getdleta(i) = EvaluateLagrangePolyDerivative( i, eta, self % Ny , self % eta)
         end do

      end function NodalStorage_getdleta

      function NodalStorage_getdlzeta(self, zeta)
         implicit none
         class(NodalStorage), intent(in)  :: self
         real(kind=RP),       intent(in)  :: zeta
         real(kind=RP)                    :: NodalStorage_getdlzeta(0:self % Nz)
!
!        ---------------
!        Local variables
!        ---------------
!
         integer  :: i

         do i = 0 , self % Nz
            NodalStorage_getdlzeta(i) = EvaluateLagrangePolyDerivative( i, zeta, self % Nz , self % zeta)
         end do

      end function NodalStorage_getdlzeta

END Module NodalStorageClass
