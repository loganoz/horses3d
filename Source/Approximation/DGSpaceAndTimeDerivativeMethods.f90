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
      MODULE DGTimeDerivativeMethods
      USE SMConstants
      use DGInviscidDiscretization
      use DGViscousDiscretization
      use DGWeakIntegrals
!
!     ========      
      CONTAINS 
!     ========      
!
!////////////////////////////////////////////////////////////////////////
!
      subroutine Initialize_SpaceAndTimeMethods()
         use PhysicsStorage
         implicit none
!
!        TODO: this should be selected from the paramfile options
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
         type(NodalStorage)         :: spA
         real(kind=RP)              :: t
!
!        ---------------
!        Local variables
!        ---------------
!
         integer     :: eID , iVar 
!
!$omp barrier
!$omp do schedule(runtime)
         do eID = 1 , size(mesh % elements)
!
!           Perform volume integrals
!           ------------------------            
            call TimeDerivative_VolumetricContribution( mesh % elements(eID) , spA , t)
!
!           Perform surface integrals
!           -------------------------
            call TimeDerivative_FacesContribution( mesh % elements(eID) , spA , t)
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
         real ( kind=RP ) :: inviscidContravariantFlux ( 0:spA % N , 0:spA % N , 0:spA % N , 1 : N_EQN , 1:NDIM ) 
         real ( kind=RP ) :: viscousContravariantFlux  ( 0:spA % N , 0:spA % N , 0:spA % N , 1 : N_EQN , 1:NDIM ) 
         real ( kind=RP ) :: contravariantFlux         ( 0:spA % N , 0:spA % N , 0:spA % N , 1 : N_EQN , 1:NDIM ) 
         integer          :: eID

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
         type(NodalStorage), intent(in) :: spA
         real(kind=RP),      intent(in) :: time
         EXTERNAL                       :: externalStateProcedure
         EXTERNAL                       :: externalGradientsProcedure

         call ViscousMethod % ComputeGradient( mesh , spA , time , externalStateProcedure , externalGradientsProcedure )

      end subroutine DGSpatial_ComputeGradient
!
!////////////////////////////////////////////////////////////////////////////////////////
!
   END MODULE DGTimeDerivativeMethods
