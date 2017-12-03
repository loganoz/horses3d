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
#include "Includes.h"
module SpatialDiscretization
      use SMConstants
      use DGInviscidDiscretization
      use DGViscousDiscretization
      use DGWeakIntegrals
      use MeshTypes
!
!     ========      
      CONTAINS 
!     ========      
!
!///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
!
      subroutine Initialize_SpaceAndTimeMethods(controlVariables)
         use PhysicsStorage
         use FTValueDictionaryClass
         use mainKeywordsModule
         use Headers
         use MPI_Process_Info
         implicit none
         class(FTValueDictionary),  intent(in)  :: controlVariables
         character(len=LINE_LENGTH)       :: inviscidDiscretization
         interface
            subroutine toLower(str)
               character(*), intent(in out) :: str
            end subroutine toLower
         end interface

         if ( MPI_Process % isRoot ) then
            write(STD_OUT,'(/)')
            call Section_Header("Spatial discretization scheme")
            write(STD_OUT,'(/)')
         end if

         inviscidDiscretization = controlVariables % stringValueForKey(inviscidDiscretizationKey,requestedLength = LINE_LENGTH)

         call toLower(inviscidDiscretization)
      
         select case ( trim(inviscidDiscretization) )

         case ( "standard" )
            if (.not. allocated(InviscidMethod)) allocate( StandardDG_t  :: InviscidMethod )

         case ( "split-form")
            if (.not. allocated(InviscidMethod)) allocate(SplitDG_t   :: InviscidMethod)

         case default
            write(STD_OUT,'(A,A,A)') 'Requested inviscid discretization "',trim(inviscidDiscretization),'" is not implemented.'
            write(STD_OUT,'(A)') "Implemented discretizations are:"
            write(STD_OUT,'(A)') "  * Standard"
            write(STD_OUT,'(A)') "  * Split-Form"
            errorMessage(STD_OUT)
            stop 

         end select
            
         call InviscidMethod % Initialize(controlVariables)
         
         if ( flowIsNavierStokes ) then
            if (.not. allocated(ViscousMethod)) allocate( BassiRebay1_t :: ViscousMethod  ) 
   
         else
            if (.not. allocated(ViscousMethod)) allocate( ViscousMethod_t  :: ViscousMethod )
            
         end if
         
      end subroutine Initialize_SpaceAndTimeMethods
!
!////////////////////////////////////////////////////////////////////////
!
      subroutine TimeDerivative_ComputeQDot( mesh , t )
         use HexMeshClass
         use ElementClass
         use PhysicsStorage
         implicit none
         type(HexMesh)              :: mesh
         real(kind=RP)              :: t
!
!        ---------------
!        Local variables
!        ---------------
!
         integer     :: eID , i, j, k
         integer     :: Nx, Ny, Nz
         interface
            subroutine UserDefinedSourceTerm(mesh, time, thermodynamics_, dimensionless_, refValues_)
               USE HexMeshClass
               use PhysicsStorage
               IMPLICIT NONE
               CLASS(HexMesh)                        :: mesh
               REAL(KIND=RP)                         :: time
               type(Thermodynamics_t),    intent(in) :: thermodynamics_
               type(Dimensionless_t),     intent(in) :: dimensionless_
               type(RefValues_t),         intent(in) :: refValues_
            end subroutine UserDefinedSourceTerm
         end interface

!$omp do schedule(runtime) 
         do eID = 1 , size(mesh % elements)
            Nx = mesh % elements(eID) % Nxyz(1)
            Ny = mesh % elements(eID) % Nxyz(2)
            Nz = mesh % elements(eID) % Nxyz(3)
!
!           Perform volume integrals
!           ------------------------            
            call TimeDerivative_VolumetricContribution( mesh % elements(eID) , t)
!
!           Perform surface integrals
!           -------------------------
            call TimeDerivative_FacesContribution( mesh % elements(eID) , t, mesh)

!
!           Scale with the Jacobian
!           -----------------------
            do k = 0, Nz   ; do j = 0, Ny    ; do i = 0, Nx
               mesh % elements(eID) % storage % QDot(:,i,j,k) = mesh % elements(eID) % storage % QDot(:,i,j,k) &
                                          / mesh % elements(eID) % geom % jacobian(i,j,k)
            end do         ; end do          ; end do
         end do
!$omp end do
!
!        Add a source term
!        -----------------
         call UserDefinedSourceTerm(mesh, t, thermodynamics, dimensionless, refValues)

      end subroutine TimeDerivative_ComputeQDot
!
!///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
!
!     --------------------------
!     TOdo: Add description here
!     --------------------------
!
      subroutine TimeDerivative_VolumetricContribution( e , t )
         use HexMeshClass
         use ElementClass
         use PhysicsStorage
         implicit none
         type(Element)      :: e
         real(kind=RP)      :: t

!
!        ---------------
!        Local variables
!        ---------------
!
         real(kind=RP) :: inviscidContravariantFlux ( 1:NCONS, 0:e%Nxyz(1) , 0:e%Nxyz(2) , 0:e%Nxyz(3), 1:NDIM ) 
         real(kind=RP) :: fSharp(1:NCONS, 0:e%Nxyz(1), 0:e%Nxyz(1), 0:e%Nxyz(2), 0:e%Nxyz(3))
         real(kind=RP) :: gSharp(1:NCONS, 0:e%Nxyz(2), 0:e%Nxyz(1), 0:e%Nxyz(2), 0:e%Nxyz(3))
         real(kind=RP) :: hSharp(1:NCONS, 0:e%Nxyz(3), 0:e%Nxyz(1), 0:e%Nxyz(2), 0:e%Nxyz(3))
         real(kind=RP) :: viscousContravariantFlux  ( 1:NCONS, 0:e%Nxyz(1) , 0:e%Nxyz(2) , 0:e%Nxyz(3), 1:NDIM ) 
         real(kind=RP) :: contravariantFlux         ( 1:NCONS, 0:e%Nxyz(1) , 0:e%Nxyz(2) , 0:e%Nxyz(3), 1:NDIM ) 
         integer       :: eID
!
!        Compute inviscid and viscous contravariant fluxes
!        -------------------------------------------------
         call InviscidMethod % ComputeInnerFluxes ( e , inviscidContravariantFlux ) 
         call ViscousMethod  % ComputeInnerFluxes ( e , viscousContravariantFlux  ) 

         select type ( InviscidMethod )
         type is (StandardDG_t)
!
!           Compute the total Navier-Stokes flux
!           ------------------------------------
            if ( flowIsNavierStokes ) then
               contravariantFlux = inviscidContravariantFlux - viscousContravariantFlux
            else
               contravariantFlux = inviscidContravariantFlux
            end if
!
!           Perform the Weak Volume Green integral
!           --------------------------------------
            e % storage % QDot = ScalarWeakIntegrals % StdVolumeGreen ( e , contravariantFlux ) 

         type is (SplitDG_t)
!
!           Compute sharp fluxes for skew-symmetric approximations
!           ------------------------------------------------------
            call InviscidMethod % ComputeSplitFormFluxes(e, inviscidContravariantFlux, fSharp, gSharp, hSharp)
!
!           Peform the Weak volume green integral
!           -------------------------------------
            if ( .not. flowIsNavierStokes ) viscousContravariantFlux = 0.0_RP
            e % storage % QDot = -ScalarWeakIntegrals % SplitVolumeDivergence( e, fSharp, gSharp, hSharp, viscousContravariantFlux)

         end select

      end subroutine TimeDerivative_VolumetricContribution
!
!///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
!
!     --------------------------
!     TOdo: Add description here
!     --------------------------
      subroutine TimeDerivative_FacesContribution( e , t , mesh)
         use HexMeshClass
         use PhysicsStorage
         implicit none
         type(Element)           :: e
         real(kind=RP)           :: t
         type(HexMesh)           :: mesh

         e % storage % QDot = e % storage % QDot - ScalarWeakIntegrals % StdFace( e, &
                      mesh % faces(e % faceIDs(EFRONT))  % storage(e % faceSide(EFRONT))  % fStar, &
                      mesh % faces(e % faceIDs(EBACK))   % storage(e % faceSide(EBACK))   % fStar, &
                      mesh % faces(e % faceIDs(EBOTTOM)) % storage(e % faceSide(EBOTTOM)) % fStar, &
                      mesh % faces(e % faceIDs(ERIGHT))  % storage(e % faceSide(ERIGHT))  % fStar, &
                      mesh % faces(e % faceIDs(ETOP))    % storage(e % faceSide(ETOP))    % fStar, &
                      mesh % faces(e % faceIDs(ELEFT))   % storage(e % faceSide(ELEFT))   % fStar )

      end subroutine TimeDerivative_FacesContribution
!
!////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
!
!              GRADIENT PROCEDURES
!              -------------------
!
!////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
!
      subroutine DGSpatial_ComputeGradient( mesh , spA, time , externalStateProcedure , externalGradientsProcedure )
         use HexMeshClass
         use PhysicsStorage
         implicit none
         type(HexMesh)                  :: mesh
         type(NodalStorage)             :: spA(0:)
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

         call ViscousMethod % ComputeGradient( mesh , spA, time , externalStateProcedure , externalGradientsProcedure )

      end subroutine DGSpatial_ComputeGradient
!
!////////////////////////////////////////////////////////////////////////////////////////
!
end module SpatialDiscretization
