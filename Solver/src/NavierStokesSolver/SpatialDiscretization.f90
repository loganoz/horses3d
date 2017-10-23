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
         integer     :: eID , i, j, k
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
            do k = 0, Nz   ; do j = 0, Ny    ; do i = 0, Nx
               mesh % elements(eID) % storage % QDot(:,i,j,k) = mesh % elements(eID) % storage % QDot(:,i,j,k) &
                                          / mesh % elements(eID) % geom % jacobian(i,j,k)
            end do         ; end do          ; end do
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
         real(kind=RP) :: inviscidContravariantFlux ( 1:NCONS, 0:spA % Nx , 0:spA % Ny , 0:spA % Nz, 1:NDIM ) 
         real(kind=RP) :: viscousContravariantFlux  ( 1:NCONS, 0:spA % Nx , 0:spA % Ny , 0:spA % Nz, 1:NDIM ) 
         real(kind=RP) :: contravariantFlux         ( 1:NCONS, 0:spA % Nx , 0:spA % Ny , 0:spA % Nz, 1:NDIM ) 
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
         e % storage % QDot = ScalarWeakIntegrals % StdVolumeGreen ( e , spA , contravariantFlux ) 

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

         e % storage % QDot = e % storage % QDot - ScalarWeakIntegrals % StdFace( e, spA, e % storage % Fstarb(:,:,:,:) ) 

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
            subroutine externalStateProcedure(x,t,nHat,Q,boundaryName)
               use SMConstants
               real(kind=RP)   , intent(in)    :: x(3), t, nHat(3)
               real(kind=RP)   , intent(inout) :: Q(:)
               character(len=*), intent(in)    :: boundaryName
            end subroutine externalStateProcedure
            
            subroutine externalGradientsProcedure(x,t,nHat,gradU,boundaryName)
               use SMConstants
               real(kind=RP)   , intent(in)    :: x(3), t, nHat(3)
               real(kind=RP)   , intent(inout) :: gradU(:,:)
               character(len=*), intent(in)    :: boundaryName
            end subroutine externalGradientsProcedure
         end interface

         call ViscousMethod % ComputeGradient( mesh , spA , time , externalStateProcedure , externalGradientsProcedure )

      end subroutine DGSpatial_ComputeGradient
!
!////////////////////////////////////////////////////////////////////////////////////////
!
end module SpatialDiscretization
