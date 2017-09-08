!
!////////////////////////////////////////////////////////////////////////
!
!      DGTimeDerivativeRoutines.f95
!      Created: 2008-07-13 16:13:12 -0400 
!      By: David Kopriva  
!
!      3D version by D.A. Kopriva 6/17/15, 12:35 PM
!
!
!////////////////////////////////////////////////////////////////////////////////////////
!
module SpatialDiscretization
      use SMConstants
      use DGInviscidDiscretization
      use DGViscousDiscretization
      use DGWeakIntegrals
!
!     ========      
      CONTAINS 
!     ========      
!
!///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
!
      subroutine Initialize_SpaceAndTimeMethods()
         use PhysicsStorage
         implicit none
!
!        TOdo: this should be selected from the paramfile options
         allocate( StandardDG_t  :: InviscidMethod )
         
         if ( flowIsNavierStokes ) then
            allocate( BassiRebay1_t :: ViscousMethod  ) 
   
         else
            allocate( ViscousMethod_t  :: ViscousMethod )
            
         end if
         
      end subroutine Initialize_SpaceAndTimeMethods
!
!////////////////////////////////////////////////////////////////////////
!
      subroutine TimeDerivative_ComputeQDot( mesh , spA , t )
         use HexMeshClass
         use ElementClass
         use NodalStorageClass
         use PhysicsStorage
         implicit none
         type(HexMesh)              :: mesh
         type(NodalStorage)         :: spA(0:,0:,0:)
         real(kind=RP)              :: t
!
!        ---------------
!        Local variables
!        ---------------
!
         integer     :: eID , iVar 
         integer     :: Nx, Ny, Nz
!
!$omp barrier
!$omp do schedule(runtime)
         do eID = 1 , size(mesh % elements)
            Nx = mesh % elements(eID) % Nxyz(1)
            Ny = mesh % elements(eID) % Nxyz(2)
            Nz = mesh % elements(eID) % Nxyz(3)
!
!           Perform volume integrals
!           ------------------------            
            call TimeDerivative_VolumetricContribution( mesh % elements(eID) , spA(Nx,Ny,Nz) , t)
!
!           Perform surface integrals
!           -------------------------
            call TimeDerivative_FacesContribution( mesh % elements(eID) , spA(Nx,Ny,Nz) , t)

!
!           Scale with the Jacobian
!           -----------------------
            do iVar = 1 , N_EQN
               mesh % elements(eID) % QDot(:,:,:,iVar) = mesh % elements(eID) % QDot(:,:,:,iVar) &
                                          / mesh % elements(eID) % geom % jacobian
            end do
         end do
!$omp end do
      end subroutine TimeDerivative_ComputeQDot
!
!///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
!
!     --------------------------
!     TOdo: Add description here
!     --------------------------
!
      subroutine TimeDerivative_VolumetricContribution( e , spA , t )
         use HexMeshClass
         use ElementClass
         use NodalStorageClass
         use PhysicsStorage
         implicit none
         type(Element)      :: e
         type(NodalStorage) :: spA
         real(kind=RP)      :: t

!
!        ---------------
!        Local variables
!        ---------------
!
         real(kind=RP) :: inviscidContravariantFlux ( 0:spA % Nx , 0:spA % Ny , 0:spA % Nz , 1 : N_EQN , 1:NDIM ) 
         real(kind=RP) :: viscousContravariantFlux  ( 0:spA % Nx , 0:spA % Ny , 0:spA % Nz , 1 : N_EQN , 1:NDIM ) 
         real(kind=RP) :: contravariantFlux         ( 0:spA % Nx , 0:spA % Ny , 0:spA % Nz , 1 : N_EQN , 1:NDIM ) 
         integer       :: eID

!
!        Compute inviscid and viscous contravariant fluxes
!        -------------------------------------------------
         call InviscidMethod % ComputeInnerFluxes ( e , spA , inviscidContravariantFlux ) 
         call ViscousMethod  % ComputeInnerFluxes ( e , spA , viscousContravariantFlux  ) 
!
!        Compute the total Navier-Stokes flux
!        ------------------------------------
         if ( flowIsNavierStokes ) then
            contravariantFlux = inviscidContravariantFlux - viscousContravariantFlux
         else
            contravariantFlux = inviscidContravariantFlux
         end if
!
!        Perform the Weak Volume Green integral
!        --------------------------------------
         e % QDot = ScalarWeakIntegrals % StdVolumeGreen ( N_EQN , e , spA , contravariantFlux ) 

      end subroutine TimeDerivative_VolumetricContribution
!
!///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
!
!     --------------------------
!     TOdo: Add description here
!     --------------------------
      subroutine TimeDerivative_FacesContribution( e , spA , t )
         use HexMeshClass
         use NodalStorageClass
         use PhysicsStorage
         implicit none
         type(Element)           :: e
         type(NodalStorage)      :: spA
         real(kind=RP)           :: t
!
!        LEFT face
!        ---------
         e % QDot = e % QDot - ScalarWeakIntegrals % StdFace( e , spA , ELEFT , e % Fstarb(:,:,:,ELEFT) ) 
!
!        RIGHT face
!        ---------
         e % QDot = e % QDot - ScalarWeakIntegrals % StdFace( e , spA , ERIGHT , e % Fstarb(:,:,:,ERIGHT) )
!
!        BOTTOM face
!        ---------
         e % QDot = e % QDot - ScalarWeakIntegrals % StdFace( e , spA , EBOTTOM , e % Fstarb(:,:,:,EBOTTOM) )
!
!        TOP face
!        ---------
         e % QDot = e % QDot - ScalarWeakIntegrals % StdFace( e , spA , ETOP , e % Fstarb(:,:,:,ETOP) )
!
!        BACK face
!        ---------
         e % QDot = e % QDot - ScalarWeakIntegrals % StdFace( e , spA , EBACK , e % Fstarb(:,:,:,EBACK) )
!
!        FRONT face
!        ---------
         e % QDot = e % QDot - ScalarWeakIntegrals % StdFace( e , spA , EFRONT , e % Fstarb(:,:,:,EFRONT) )

      end subroutine TimeDerivative_FacesContribution
!
!////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
!
!              GRADIENT PROCEDURES
!              -------------------
!
!////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
!
      subroutine DGSpatial_ComputeGradient( mesh , spA , time , externalStateProcedure , externalGradientsProcedure )
         use HexMeshClass
         use NodalStorageClass
         use PhysicsStorage
         implicit none
         type(HexMesh)                  :: mesh
         type(NodalStorage), intent(in) :: spA(0:,0:,0:)
         real(kind=RP),      intent(in) :: time
         interface
            subroutine externalStateSubroutine(x,t,nHat,Q,boundaryName)
               use SMConstants
               real(kind=RP)   , intent(in)    :: x(3), t, nHat(3)
               real(kind=RP)   , intent(inout) :: Q(:)
               character(len=*), intent(in)    :: boundaryName
            end subroutine externalStateSubroutine
            
            subroutine externalGradientsSubroutine(x,t,nHat,gradU,boundaryName)
               use SMConstants
               real(kind=RP)   , intent(in)    :: x(3), t, nHat(3)
               real(kind=RP)   , intent(inout) :: gradU(:,:)
               character(len=*), intent(in)    :: boundaryName
            end subroutine externalGradientsSubroutine
         end interface

         call ViscousMethod % ComputeGradient( mesh , spA , time , externalStateProcedure , externalGradientsProcedure )

      end subroutine DGSpatial_ComputeGradient
!
!////////////////////////////////////////////////////////////////////////////////////////
!
   !! Old routine... TOdo: See what happens with this code
!~   subroutine DGGradSpaceDerivative( u, uL, uR,  N, D, b, deriv ) 
!~!
!~!  ----------------------------------------------------------
!~!  The one dimensional space derivative for the DG method for
!~!  a vector state. This is Algorithm 92 in the book.
!~!  ----------------------------------------------------------
!~!
!~      use PhysicsStorage
!~      use PolynomialInterpAndDerivsModule
!~      IMPLICIT NONE
!~!
!~!     ---------
!~!     Arguments
!~!     ---------
!~!
!~      INTEGER                                   :: N
!~      REAL(KIND=RP), DIMENSION(0:N, N_GRAD_EQN) :: u
!~      REAL(KIND=RP), DIMENSION(N_GRAD_EQN)      :: uL, uR
!~      REAL(KIND=RP), DIMENSION(0:N, 0:N)        :: D
!~      REAL(KIND=RP), DIMENSION(0:N, 2)          :: b
!~!
!~!     -----------------
!~!     Output variables:
!~!     -----------------
!~!
!~      REAL(KIND=RP), DIMENSION(0:N,N_GRAD_EQN), INTENT(OUT) :: deriv
!~!
!~!     ---------------
!~!     Local variables
!~!     ---------------
!~!
!~      INTEGER            :: j, k
!~      INTEGER, PARAMETER :: LEFT = 1, RIGHT = 2
!~!
!~!     ----------------------------
!~!     Internal points contribution
!~!     ----------------------------
!~!
!~      do k = 1, N_GRAD_EQN
!~         CALL PolyDirectMatrixMultiplyDeriv( u(:,k), deriv(:,k), D, N )
!~      end do
!~!
!~!     ----------------------------
!~!     Boundary points contribution
!~!     ----------------------------
!~!
!~      do j = 0,N  
!~         deriv(j,:) = deriv(j,:) + uR*b(j,RIGHT) + uL*b(j,LEFT)
!~      end do
     
!~   end subroutine DGGradSpaceDerivative   
!
!
!///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
!!
!      subroutine ComputeContravariantFlux( e, contravariantFlux )
!!
!!     --------------------------------------
!!     As described, compute
!      
!!     \[
!!        \tilde f^i = J\vec a^i \cdot \vec F
!!     \]
!!     --------------------------------------
!!      
!      use ElementClass
!      use PhysicsStorage
!      use Physics
!      IMPLICIT NONE
!!
!!     -----------------
!!     Input parameters:
!!     -----------------
!!
!      TYPE(Element)                          :: e
!      REAL(KIND=RP), dimension(0:,0:,0:,:,:) :: contravariantFlux 
!!
!!     ---------------
!!     Local variables
!!     ---------------
!!
!      INTEGER                         :: n, m, l, nv
!      REAL(KIND=RP), DIMENSION(N_EQN) :: ff, gg, hh
!      
!      do l = 0, e % Nxyz(3)
!         do m = 0, e % Nxyz(2)
!            do n = 0, e % Nxyz(1)
!               do nv = 1, N_EQN
!                  contravariantFlux(n,m,l,nv,1) = e % geom % jGradXi(1,n,m,l)  *ff(nv) +   &
!                                                  e % geom % jGradXi(2,n,m,l)  *gg(nv) +   &
!                                                  e % geom % jGradXi(3,n,m,l)  *hh(nv)
!                  contravariantFlux(n,m,l,nv,2) = e % geom % jGradEta(1,n,m,l) *ff(nv) +   &
!                                                  e % geom % jGradEta(2,n,m,l) *gg(nv) +   &
!                                                  e % geom % jGradEta(3,n,m,l) *hh(nv)
!                  contravariantFlux(n,m,l,nv,3) = e % geom % jGradZeta(1,n,m,l)*ff(nv) +   &
!                                                  e % geom % jGradZeta(2,n,m,l)*gg(nv) +   &
!                                                  e % geom % jGradZeta(3,n,m,l)*hh(nv)
!               end do
!               
!            end do
!         end do
!      end do
!    
!      end subroutine ComputeContravariantFlux
!
!////////////////////////////////////////////////////////////////////////
!
end module SpatialDiscretization
