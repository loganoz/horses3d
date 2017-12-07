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
#ifdef _HAS_MPI_
      use mpi
#endif
!
!     ========      
      CONTAINS 
!     ========      
!
!////////////////////////////////////////////////////////////////////////////////////////
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
      subroutine TimeDerivative_ComputeQDot( mesh , t, externalState, externalGradients )
         use HexMeshClass
         use ElementClass
         use PhysicsStorage
         use MPI_Face_Class
         use MPI_Process_Info
         implicit none
         type(HexMesh)              :: mesh
         real(kind=RP)              :: t
         external                   :: externalState, externalGradients
!
!        ---------------
!        Local variables
!        ---------------
!
         integer     :: eID , i, j, k, ierr, fID
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
!
!        ****************
!        Volume integrals
!        ****************
!
!$omp do schedule(runtime) 
         do eID = 1 , size(mesh % elements)
            call TimeDerivative_VolumetricContribution( mesh % elements(eID) , t)
         end do
!$omp end do nowait
!
!        ******************************************
!        Compute Riemann solver of non-shared faces
!        ******************************************
!
!$omp do schedule(runtime) 
         do fID = 1, size(mesh % faces) 
            associate( f => mesh % faces(fID)) 
            select case (f % faceType) 
            case (HMESH_INTERIOR) 
               CALL computeElementInterfaceFlux ( f ) 
 
            case (HMESH_BOUNDARY) 
               CALL computeBoundaryFlux(f, t, externalState, externalGradients) 
 
            end select 
            end associate 
         end do 
!$omp end do 
!
!        ***************************************************************
!        Surface integrals and scaling of elements with non-shared faces
!        ***************************************************************
! 
!$omp do schedule(runtime) 
         do eID = 1, size(mesh % elements) 
            associate(e => mesh % elements(eID)) 
            if ( e % hasSharedFaces ) cycle
            call TimeDerivative_FacesContribution(e, t, mesh) 
 
            do k = 0, e % Nxyz(3) ; do j = 0, e % Nxyz(2) ; do i = 0, e % Nxyz(1) 
               e % storage % QDot(:,i,j,k) = e % storage % QDot(:,i,j,k) / e % geom % jacobian(i,j,k) 
            end do         ; end do          ; end do 
            end associate 
         end do
!$omp end do
!
!        ****************************
!        Wait until messages are sent
!        ****************************
!
!$omp single
         if ( MPI_Process % doMPIAction ) then
            if ( flowIsNavierStokes ) then 
               call WaitUntilGradientsAreReady(mpi_faces) 
            else  
               call WaitUntilSolutionIsReady(mpi_faces) 
            end if          
         end if
!$omp end single
!
!        **************************************
!        Compute Riemann solver of shared faces
!        **************************************
!
!$omp do schedule(runtime) 
         do fID = 1, size(mesh % faces) 
            associate( f => mesh % faces(fID)) 
            select case (f % faceType) 
            case (HMESH_MPI) 
               CALL computeMPIFaceFlux ( f ) 
            end select 
            end associate 
         end do 
!$omp end do 
!
!        ***********************************************************
!        Surface integrals and scaling of elements with shared faces
!        ***********************************************************
! 
!$omp do schedule(runtime) 
         do eID = 1, size(mesh % elements) 
            associate(e => mesh % elements(eID)) 
            if ( .not. e % hasSharedFaces ) cycle
            call TimeDerivative_FacesContribution(e, t, mesh) 
 
            do k = 0, e % Nxyz(3) ; do j = 0, e % Nxyz(2) ; do i = 0, e % Nxyz(1) 
               e % storage % QDot(:,i,j,k) = e % storage % QDot(:,i,j,k) / e % geom % jacobian(i,j,k) 
            end do         ; end do          ; end do 
            end associate 
         end do
!$omp end do
!
!        Add a source term
!        -----------------
         call UserDefinedSourceTerm(mesh, t, thermodynamics, dimensionless, refValues)
!
!        Add a MPI Barrier
!        -----------------
#ifdef _HAS_MPI_
!$omp single
         if ( MPI_Process % doMPIAction ) call mpi_barrier(MPI_COMM_WORLD, ierr)
!$omp end single
#endif

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
!///////////////////////////////////////////////////////////////////////////////////////////// 
! 
!        Riemann solver drivers 
!        ---------------------- 
! 
!///////////////////////////////////////////////////////////////////////////////////////////// 
! 
      SUBROUTINE computeElementInterfaceFlux(f)
         use FaceClass
         use Physics
         use PhysicsStorage
         IMPLICIT NONE
         TYPE(Face)   , INTENT(inout) :: f   
         integer       :: i, j
         real(kind=RP) :: inv_flux(1:N_EQN,0:f % Nf(1),0:f % Nf(2))
         real(kind=RP) :: visc_flux(1:N_EQN,0:f % Nf(1),0:f % Nf(2))
         real(kind=RP) :: flux(1:N_EQN,0:f % Nf(1),0:f % Nf(2))
!         
!
!        --------------
!        Invscid fluxes:
!        --------------
!
         DO j = 0, f % Nf(2)
            DO i = 0, f % Nf(1)
               CALL RiemannSolver(QLeft  = f % storage(1) % Q(:,i,j), &
                                  QRight = f % storage(2) % Q(:,i,j), &
                                  nHat   = f % geom % normal(:,i,j), &
                                  flux   = inv_flux(:,i,j) )

               CALL ViscousMethod % RiemannSolver(QLeft = f % storage(1) % Q(:,i,j), &
                                                  QRight = f % storage(2) % Q(:,i,j), &
                                                  U_xLeft = f % storage(1) % U_x(:,i,j), &
                                                  U_yLeft = f % storage(1) % U_y(:,i,j), &
                                                  U_zLeft = f % storage(1) % U_z(:,i,j), &
                                                  U_xRight = f % storage(2) % U_x(:,i,j), &
                                                  U_yRight = f % storage(2) % U_y(:,i,j), &
                                                  U_zRight = f % storage(2) % U_z(:,i,j), &
                                                  nHat = f % geom % normal(:,i,j) , &
                                                  flux  = visc_flux(:,i,j) )
               
!
!              Multiply by the Jacobian
!              ------------------------
               flux(:,i,j) = ( inv_flux(:,i,j) - visc_flux(:,i,j) ) * f % geom % jacobian(i,j)
               
            END DO   
         END DO  
!
!        ---------------------------
!        Return the flux to elements
!        ---------------------------
!
         call f % ProjectFluxToElements(flux, (/1,2/))

      END SUBROUTINE computeElementInterfaceFlux

      SUBROUTINE computeMPIFaceFlux(f)
         use FaceClass
         use Physics
         use PhysicsStorage
         IMPLICIT NONE
         TYPE(Face)   , INTENT(inout) :: f   
         integer       :: i, j
         integer       :: thisSide
         real(kind=RP) :: inv_flux(1:N_EQN,0:f % Nf(1),0:f % Nf(2))
         real(kind=RP) :: visc_flux(1:N_EQN,0:f % Nf(1),0:f % Nf(2))
         real(kind=RP) :: flux(1:N_EQN,0:f % Nf(1),0:f % Nf(2))
!         
!
!        --------------
!        Invscid fluxes:
!        --------------
!
         DO j = 0, f % Nf(2)
            DO i = 0, f % Nf(1)
               CALL RiemannSolver(QLeft  = f % storage(1) % Q(:,i,j), &
                                  QRight = f % storage(2) % Q(:,i,j), &
                                  nHat   = f % geom % normal(:,i,j), &
                                  flux   = inv_flux(:,i,j) )

               CALL ViscousMethod % RiemannSolver(QLeft = f % storage(1) % Q(:,i,j), &
                                                  QRight = f % storage(2) % Q(:,i,j), &
                                                  U_xLeft = f % storage(1) % U_x(:,i,j), &
                                                  U_yLeft = f % storage(1) % U_y(:,i,j), &
                                                  U_zLeft = f % storage(1) % U_z(:,i,j), &
                                                  U_xRight = f % storage(2) % U_x(:,i,j), &
                                                  U_yRight = f % storage(2) % U_y(:,i,j), &
                                                  U_zRight = f % storage(2) % U_z(:,i,j), &
                                                  nHat = f % geom % normal(:,i,j) , &
                                                  flux  = visc_flux(:,i,j) )
               
!
!              Multiply by the Jacobian
!              ------------------------
               flux(:,i,j) = ( inv_flux(:,i,j) - visc_flux(:,i,j) ) * f % geom % jacobian(i,j)
               
            END DO   
         END DO  
!
!        ---------------------------
!        Return the flux to elements: The sign in eR % storage % FstarB has already been accouted.
!        ---------------------------
!
         thisSide = maxloc(f % elementIDs, dim = 1)
         call f % ProjectFluxToElements(flux, (/thisSide, HMESH_NONE/))


      end subroutine ComputeMPIFaceFlux

      SUBROUTINE computeBoundaryFlux(f, time, externalStateProcedure , externalGradientsProcedure)
      USE ElementClass
      use FaceClass
      USE DGViscousDiscretization
      USE Physics
      use PhysicsStorage
      USE BoundaryConditionFunctions
      IMPLICIT NONE
!
!     ---------
!     Arguments
!     ---------
!
      type(Face),    intent(inout) :: f
      REAL(KIND=RP)                :: time
      EXTERNAL                     :: externalStateProcedure
      EXTERNAL                     :: externalGradientsProcedure
!
!     ---------------
!     Local variables
!     ---------------
!
      INTEGER                         :: i, j
      INTEGER, DIMENSION(2)           :: N
      REAL(KIND=RP)                   :: inv_flux(N_EQN)
      REAL(KIND=RP)                   :: UGradExt(NDIM , N_GRAD_EQN) , visc_flux(N_EQN)
      real(kind=RP)                   :: fStar(N_EQN, 0:f % Nf(1), 0: f % Nf(2))
      CHARACTER(LEN=BC_STRING_LENGTH) :: boundaryType
            
      boundaryType = f % boundaryType
      
      DO j = 0, f % Nf(2)
         DO i = 0, f % Nf(1)
!
!           Inviscid part
!           -------------
            f % storage(2) % Q(:,i,j) = f % storage(1) % Q(:,i,j)
            CALL externalStateProcedure( f % geom % x(:,i,j), &
                                         time, &
                                         f % geom % normal(:,i,j), &
                                         f % storage(2) % Q(:,i,j),&
                                         boundaryType )
            CALL RiemannSolver(QLeft  = f % storage(1) % Q(:,i,j), &
                               QRight = f % storage(2) % Q(:,i,j), &   
                               nHat   = f % geom % normal(:,i,j), &
                               flux   = inv_flux)
!
!           ViscousPart
!           -----------
            if ( flowIsNavierStokes ) then

            
            UGradExt(IX,:) = f % storage(1) % U_x(:,i,j)
            UGradExt(IY,:) = f % storage(1) % U_y(:,i,j)
            UGradExt(IZ,:) = f % storage(1) % U_z(:,i,j)

            CALL externalGradientsProcedure(  f % geom % x(:,i,j), &
                                              time, &
                                              f % geom % normal(:,i,j), &
                                              UGradExt,&
                                              boundaryType )

            CALL ViscousMethod % RiemannSolver( QLeft = f % storage(1) % Q(:,i,j), &
                                                QRight = f % storage(2) % Q(:,i,j), &
                                                U_xLeft = f % storage(1) % U_x(:,i,j), &
                                                U_yLeft = f % storage(1) % U_y(:,i,j), &
                                                U_zLeft = f % storage(1) % U_z(:,i,j), &
                                                U_xRight = UGradExt(IX,:) , &
                                                U_yRight = UGradExt(IY,:) , &
                                                U_zRight = UGradExt(IZ,:) , &
                                                nHat = f % geom % normal(:,i,j), &
                                                flux = visc_flux )
            else
               visc_flux = 0.0_RP
            end if

            fStar(:,i,j) = (inv_flux - visc_flux) * f % geom % jacobian(i,j)
         END DO   
      END DO   

      call f % ProjectFluxToElements(fStar, (/1, HMESH_NONE/))

      END SUBROUTINE computeBoundaryFlux
!
!//////////////////////////////////////////////////////////////////////// 
! 

!
!////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
!
!              GRADIENT PROCEDURES
!              -------------------
!
!////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
!
      subroutine DGSpatial_ComputeGradient( mesh , time , externalStateProcedure , externalGradientsProcedure )
         use HexMeshClass
         use PhysicsStorage
         implicit none
         type(HexMesh)                  :: mesh
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

         call ViscousMethod % ComputeGradient( mesh , time , externalStateProcedure , externalGradientsProcedure )

      end subroutine DGSpatial_ComputeGradient
!
!////////////////////////////////////////////////////////////////////////////////////////
!
end module SpatialDiscretization
