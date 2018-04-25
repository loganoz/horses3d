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
      use HyperbolicDiscretizations
      use EllipticDiscretizations
      use LESModels
      use SpectralVanishingViscosity
      use DGWeakIntegrals
      use MeshTypes
      use HexMeshClass
      use ElementClass
      use PhysicsStorage
      use Physics
      use MPI_Face_Class
      use MPI_Process_Info
      use DGSEMClass
      use ParticlesClass
      use FluidData
      use VariableConversion, only: NSGradientValuesForQ_0D, NSGradientValuesForQ_3D
#ifdef _HAS_MPI_
      use mpi
#endif


      abstract interface
         SUBROUTINE computeElementInterfaceFluxF(f)
            use FaceClass
            IMPLICIT NONE
            TYPE(Face)   , INTENT(inout) :: f   
         end subroutine computeElementInterfaceFluxF

         SUBROUTINE computeMPIFaceFluxF(f)
            use FaceClass
            IMPLICIT NONE
            TYPE(Face)   , INTENT(inout) :: f   
         end subroutine computeMPIFaceFluxF

         SUBROUTINE computeBoundaryFluxF(f, time, externalStateProcedure , externalGradientsProcedure)
            use SMConstants
            use FaceClass
            IMPLICIT NONE
            type(Face),    intent(inout) :: f
            REAL(KIND=RP)                :: time
            EXTERNAL                     :: externalStateProcedure
            EXTERNAL                     :: externalGradientsProcedure
         end subroutine computeBoundaryFluxF
      end interface

      procedure(computeElementInterfaceFluxF), pointer :: computeElementInterfaceFlux => computeElementInterfaceFlux_NS
      procedure(computeMPIFaceFluxF),          pointer :: computeMPIFaceFlux          => computeMPIFaceFlux_NS
      procedure(computeBoundaryFluxF),         pointer :: computeBoundaryFlux         => computeBoundaryFlux_NS

!
!     ========      
      CONTAINS 
!     ========      
!
!////////////////////////////////////////////////////////////////////////////////////////
!
      subroutine Initialize_SpaceAndTimeMethods(controlVariables, mesh)
         use FTValueDictionaryClass
         use Utilities, only: toLower
         use mainKeywordsModule
         use Headers
         use MPI_Process_Info
         implicit none
         class(FTValueDictionary),  intent(in)  :: controlVariables
         class(HexMesh)                         :: mesh
         character(len=LINE_LENGTH)       :: inviscidDiscretization
         character(len=LINE_LENGTH)       :: viscousDiscretization
         
         if (.not. mesh % child) then ! If this is a child mesh, all these constructs were already initialized for the parent mesh
         
            if ( MPI_Process % isRoot ) then
               write(STD_OUT,'(/)')
               call Section_Header("Spatial discretization scheme")
               write(STD_OUT,'(/)')
            end if
   !
   !        Initialize inviscid discretization
   !        ----------------------------------
            inviscidDiscretization = controlVariables % stringValueForKey(inviscidDiscretizationKey,requestedLength = LINE_LENGTH)

            call toLower(inviscidDiscretization)
         
            select case ( trim(inviscidDiscretization) )

            case ( "standard" )
               if (.not. allocated(HyperbolicDiscretization)) allocate( StandardDG_t  :: HyperbolicDiscretization )

            case ( "split-form")
               if (.not. allocated(HyperbolicDiscretization)) allocate( SplitDG_t     :: HyperbolicDiscretization)

            case default
               write(STD_OUT,'(A,A,A)') 'Requested inviscid discretization "',trim(inviscidDiscretization),'" is not implemented.'
               write(STD_OUT,'(A)') "Implemented discretizations are:"
               write(STD_OUT,'(A)') "  * Standard"
               write(STD_OUT,'(A)') "  * Split-Form"
               errorMessage(STD_OUT)
               stop 

            end select
               
            call HyperbolicDiscretization % Initialize(controlVariables)
   !
   !        Initialize viscous discretization
   !        ---------------------------------         
            if ( flowIsNavierStokes ) then
               call BassiRebay1     % Initialize(controlVariables)
               call BassiRebay2     % Initialize(controlVariables)
               call InteriorPenalty % Initialize(controlVariables)

               viscousDiscretization = controlVariables % stringValueForKey(viscousDiscretizationKey, requestedLength = LINE_LENGTH)
               call toLower(viscousDiscretization)
               
               select case ( trim(viscousDiscretization) )
               case("br1")
                  EllipticDiscretization => BassiRebay1

               case("br2")
                  EllipticDiscretization => BassiRebay2

               case("ip")
                  EllipticDiscretization => InteriorPenalty

               case default
                  write(STD_OUT,'(A,A,A)') 'Requested viscous discretization "',trim(viscousDiscretization),'" is not implemented.'
                  write(STD_OUT,'(A)') "Implemented discretizations are:"
                  write(STD_OUT,'(A)') "  * BR1"
                  write(STD_OUT,'(A)') "  * BR2"
                  write(STD_OUT,'(A)') "  * IP"
                  errorMessage(STD_OUT)
                  stop 

               end select

               call EllipticDiscretization % Describe
      
            else
               if (.not. associated(EllipticDiscretization)) allocate( EllipticDiscretization_t  :: EllipticDiscretization )
               call EllipticDiscretization % Initialize(controlVariables)
               
            end if

   !
   !        Initialize models
   !        -----------------
            call InitializeLESModel(LESModel, controlVariables)
         
         end if
!
!        Compute wall distances
!        ----------------------
         call mesh % ComputeWallDistances
!
!        Initialize SVV
!        --------------
         call InitializeSVV(SVV, controlVariables, mesh)
         
         if (.not. mesh % child) then
            if ( SVV % enabled ) then
               computeElementInterfaceFlux => computeElementInterfaceFlux_SVV
               computeMPIFaceFlux          => computeMPIFaceFlux_SVV
               computeBoundaryFlux         => computeBoundaryFlux_SVV

            end if
         end if
         
      end subroutine Initialize_SpaceAndTimeMethods
!
!////////////////////////////////////////////////////////////////////////
!
      subroutine Finalize_SpaceAndTimeMethods
         implicit none
         IF ( ALLOCATED(HyperbolicDiscretization) ) DEALLOCATE( HyperbolicDiscretization )
         IF ( ALLOCATED(LESModel) )       DEALLOCATE( LESModel )
      end subroutine Finalize_SpaceAndTimeMethods
!
!////////////////////////////////////////////////////////////////////////
!
      SUBROUTINE ComputeTimeDerivative( mesh, particles, time, BCFunctions)
         IMPLICIT NONE 
!
!        ---------
!        Arguments
!        ---------
!
         TYPE(HexMesh), target           :: mesh
         type(Particles_t)               :: particles
         REAL(KIND=RP)                   :: time
         type(BCFunctions_t), intent(in) :: BCFunctions(no_of_BCsets)
!
!        ---------------
!        Local variables
!        ---------------
!
         INTEGER :: k
!
!        -----------------------------------------
!        Prolongation of the solution to the faces
!        -----------------------------------------
!
!$omp parallel shared(mesh, time)
         call mesh % ProlongSolutionToFaces(NCONS)

!        ----------------
!        Update MPI Faces
!        ----------------
!
!$omp single
         call mesh % UpdateMPIFacesSolution
!$omp end single
!
!        -----------------
!        Compute gradients
!        -----------------
!
         if ( computeGradients ) then
            CALL DGSpatial_ComputeGradient(mesh , time , BCFunctions(1) % externalState)
         end if

!$omp single
         if ( flowIsNavierStokes ) then
            call mesh % UpdateMPIFacesGradients
         end if
!$omp end single
!
!        -----------------------
!        Compute time derivative
!        -----------------------
!
         call TimeDerivative_ComputeQDot(mesh = mesh , &
                                         particles = particles, &
                                         t    = time, &
                                         externalState     = BCFunctions(1) % externalState, &
                                         externalGradients = BCFunctions(1) % externalGradients )
!$omp end parallel
!
      END SUBROUTINE ComputeTimeDerivative
!
!////////////////////////////////////////////////////////////////////////
!
!     This routine computes the time derivative element by element, without considering the Riemann Solvers
!     This is useful for estimating the isolated truncation error
!
      SUBROUTINE ComputeTimeDerivativeIsolated( mesh, particles, time, BCFunctions)
         use EllipticDiscretizationClass
         IMPLICIT NONE 
!
!        ---------
!        Arguments
!        ---------
!
         TYPE(HexMesh), target      :: mesh
         type(Particles_t)          :: particles
         REAL(KIND=RP)              :: time
         type(BCFunctions_t), intent(in)  :: BCFunctions(no_of_BCsets)
!
!        ---------------
!        Local variables
!        ---------------
!
         INTEGER :: k
!
!        -----------------------------------------
!        Prolongation of the solution to the faces
!        -----------------------------------------
!
!$omp parallel shared(mesh, time)
         call mesh % ProlongSolutionToFaces(NCONS)
!
!        -----------------------------------------------------
!        Compute LOCAL gradients and prolong them to the faces
!        -----------------------------------------------------
!
         if ( computeGradients ) then
            CALL BaseClass_ComputeGradient( EllipticDiscretization, NCONS, NGRAD, mesh , time , BCFunctions(1) % externalState, NSGradientValuesForQ_0D, NSGradientValuesForQ_3D )
!
!           The prolongation is usually done in the viscous methods, but not in the BaseClass
!           ---------------------------------------------------------------------------------
            call mesh % ProlongGradientsToFaces(NGRAD)
         end if

!
!        -----------------------
!        Compute time derivative
!        -----------------------
!
         call TimeDerivative_ComputeQDotIsolated(mesh = mesh , &
                                                 t    = time )
!$omp end parallel
!
      END SUBROUTINE ComputeTimeDerivativeIsolated

      subroutine TimeDerivative_ComputeQDot( mesh , particles, t, externalState, externalGradients )
         implicit none
         type(HexMesh)              :: mesh
         type(Particles_t)          :: particles
         real(kind=RP)              :: t
         procedure(BCState_FCN)     :: externalState
         procedure(BCGradients_FCN) :: externalGradients
!
!        ---------------
!        Local variables
!        ---------------
!
         integer     :: eID , i, j, k, ierr, fID
         interface
            subroutine UserDefinedSourceTerm(mesh, time, thermodynamics_, dimensionless_, refValues_)
               use SMConstants
               USE HexMeshClass
               use PhysicsStorage
               use FluidData
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
               CALL computeElementInterfaceFlux( f ) 
 
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
!$omp do schedule(runtime) private(i,j,k)
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
#ifdef _HAS_MPI_
         if ( MPI_Process % doMPIAction ) then
!$omp single
            if ( flowIsNavierStokes ) then 
               call mesh % GatherMPIFacesGradients
            else  
               call mesh % GatherMPIFacesSolution
            end if          
!$omp end single
!
!           **************************************
!           Compute Riemann solver of shared faces
!           **************************************
!
!$omp do schedule(runtime) 
            do fID = 1, size(mesh % faces) 
               associate( f => mesh % faces(fID)) 
               select case (f % faceType) 
               case (HMESH_MPI) 
                  CALL computeMPIFaceFlux( f ) 
               end select 
               end associate 
            end do 
!$omp end do 
!
!           ***********************************************************
!           Surface integrals and scaling of elements with shared faces
!           ***********************************************************
! 
!$omp do schedule(runtime) private(i,j,k)
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
!           Add a MPI Barrier
!           -----------------
!$omp single
            call mpi_barrier(MPI_COMM_WORLD, ierr)
!$omp end single
         end if
#endif
!
!        Add a source term
!        -----------------
         if (.not. mesh % child) call UserDefinedSourceTerm(mesh, t, thermodynamics, dimensionless, refValues)

         if (.not. mesh % child) then
            if ( particles % active ) then             
!$omp do schedule(runtime)
               do eID = 1, size(mesh % elements)
                  call particles % AddSource(mesh % elements(eID), t, thermodynamics, dimensionless, refValues)
               end do
!$omp end do
            endif 
         end if

!$omp do schedule(runtime) private(i,j,k)
         do eID = 1, mesh % no_of_elements
            associate ( e => mesh % elements(eID) )
            do k = 0, e % Nxyz(3)   ; do j = 0, e % Nxyz(2) ; do i = 0, e % Nxyz(1)
               e % storage % QDot(:,i,j,k) = e % storage % QDot(:,i,j,k) + e % storage % S(:,i,j,k)
            end do                  ; end do                ; end do
            end associate
         end do
!$omp end do
      end subroutine TimeDerivative_ComputeQDot
!
!////////////////////////////////////////////////////////////////////////
!
!     -------------------------------------------------------------------------------
!     This routine computes Qdot neglecting the interaction with neighboring elements
!     and boundaries. Therefore, the external states are not needed.
!     -------------------------------------------------------------------------------
      subroutine TimeDerivative_ComputeQDotIsolated( mesh , t )
         implicit none
         type(HexMesh)              :: mesh
         real(kind=RP)              :: t
!
!        ---------------
!        Local variables
!        ---------------
!
         integer     :: eID , i, j, k, fID
         interface
            subroutine UserDefinedSourceTerm(mesh, time, thermodynamics_, dimensionless_, refValues_)
               use SMConstants
               USE HexMeshClass
               use PhysicsStorage
               use FluidData
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
!$omp end do
!
!        ******************************
!        Compute isolated flux on faces
!           Only NS version... TODO: implement SVV version
!        ******************************
!
!$omp do schedule(runtime) 
         do fID = 1, size(mesh % faces) 
            associate( f => mesh % faces(fID))
            select case (f % faceType) 
               case (HMESH_INTERIOR) 
                  call computeIsolatedFaceFluxes_NS( f )
               case (HMESH_BOUNDARY)
                  call computeIsolatedFaceFlux_NS  ( f , 1 )
               case (HMESH_MPI)
                  call computeIsolatedFaceFlux_NS  ( f , maxloc(f % elementIDs, dim = 1) )
            end select
            end associate 
         end do 
!$omp end do 
!
!        ***********************************************************
!        Surface integrals (with local data) and scaling of elements
!        ***********************************************************
! 
!$omp do schedule(runtime) private(i,j,k)
         do eID = 1, size(mesh % elements) 
            associate(e => mesh % elements(eID))
            call TimeDerivative_FacesContribution(e, t, mesh) 

            do k = 0, e % Nxyz(3) ; do j = 0, e % Nxyz(2) ; do i = 0, e % Nxyz(1) 
               e % storage % QDot(:,i,j,k) = e % storage % QDot(:,i,j,k) / e % geom % jacobian(i,j,k) 
            end do         ; end do          ; end do 
            end associate 
         end do
!$omp end do
!
!        *****************
!        Add a source term
!        *****************
!
         if (.not. mesh % child) call UserDefinedSourceTerm(mesh, t, thermodynamics, dimensionless, refValues)
!$omp do schedule(runtime) private(i,j,k)
         do eID = 1, mesh % no_of_elements
            associate ( e => mesh % elements(eID) )
            do k = 0, e % Nxyz(3)   ; do j = 0, e % Nxyz(2) ; do i = 0, e % Nxyz(1)
               e % storage % QDot(:,i,j,k) = e % storage % QDot(:,i,j,k) + e % storage % S(:,i,j,k)
            end do                  ; end do                ; end do
            end associate
         end do
!$omp end do
         
      end subroutine TimeDerivative_ComputeQDotIsolated
!
!///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
!
      subroutine TimeDerivative_VolumetricContribution( e , t )
         use HexMeshClass
         use ElementClass
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
         real(kind=RP) :: SVVContravariantFlux  ( 1:NCONS, 0:e%Nxyz(1) , 0:e%Nxyz(2) , 0:e%Nxyz(3), 1:NDIM ) 
         real(kind=RP) :: contravariantFlux         ( 1:NCONS, 0:e%Nxyz(1) , 0:e%Nxyz(2) , 0:e%Nxyz(3), 1:NDIM ) 
         integer       :: eID
!
!        *************************************
!        Compute interior contravariant fluxes
!        *************************************
!
!        Compute inviscid contravariant flux
!        -----------------------------------
         call HyperbolicDiscretization % ComputeInnerFluxes ( e , EulerFlux3D, inviscidContravariantFlux ) 
!
!        Compute viscous contravariant flux
!        ----------------------------------
         if ( .not. LESModel % active ) then
!
!           Without LES model
!           -----------------
            call EllipticDiscretization  % ComputeInnerFluxes ( NCONS, NGRAD, e , ViscousFlux3D, viscousContravariantFlux) 

         else
!
!           With LES model
!           --------------
            call EllipticDiscretization  % ComputeInnerFluxesWithSGS ( e , viscousContravariantFlux  ) 

         end if
!
!        Compute the SVV dissipation
!        ---------------------------
         if ( .not. SVV % enabled ) then
            SVVcontravariantFlux = 0.0_RP
         else
            call SVV % ComputeInnerFluxes(e, ViscousFlux3D, SVVContravariantFlux)
         end if
!
!        ************************
!        Perform volume integrals
!        ************************
!
         select type ( HyperbolicDiscretization )
         type is (StandardDG_t)
!
!           Compute the total Navier-Stokes flux
!           ------------------------------------
            contravariantFlux = inviscidContravariantFlux - viscousContravariantFlux - SVVContravariantFlux
!
!           Perform the Weak Volume Green integral
!           --------------------------------------
            e % storage % QDot = ScalarWeakIntegrals % StdVolumeGreen ( e, NCONS, contravariantFlux ) 

         type is (SplitDG_t)
!
!           Compute sharp fluxes for skew-symmetric approximations
!           ------------------------------------------------------
            call HyperbolicDiscretization % ComputeSplitFormFluxes(e, inviscidContravariantFlux, fSharp, gSharp, hSharp)
!
!           Peform the Weak volume green integral
!           -------------------------------------
            viscousContravariantFlux = viscousContravariantFlux + SVVContravariantFlux

            e % storage % QDot = -ScalarWeakIntegrals % SplitVolumeDivergence( e, fSharp, gSharp, hSharp, viscousContravariantFlux)

         end select

      end subroutine TimeDerivative_VolumetricContribution
!
!///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
!
      subroutine TimeDerivative_FacesContribution( e , t , mesh)
         use HexMeshClass
         implicit none
         type(Element)           :: e
         real(kind=RP)           :: t
         type(HexMesh)           :: mesh

         e % storage % QDot = e % storage % QDot - ScalarWeakIntegrals % StdFace( e, NCONS, &
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
      SUBROUTINE computeElementInterfaceFlux_NS(f)
         use FaceClass
         use RiemannSolvers_NS
         IMPLICIT NONE
         TYPE(Face)   , INTENT(inout) :: f   
         integer       :: i, j
         real(kind=RP) :: inv_flux(1:NCONS,0:f % Nf(1),0:f % Nf(2))
         real(kind=RP) :: visc_flux(1:NCONS,0:f % Nf(1),0:f % Nf(2))
         real(kind=RP) :: flux(1:NCONS,0:f % Nf(1),0:f % Nf(2))

         if ( .not. LESModel % active ) then
         DO j = 0, f % Nf(2)
            DO i = 0, f % Nf(1)

!      
!              --------------
!              Viscous fluxes
!              --------------
!      
               CALL EllipticDiscretization % RiemannSolver(nEqn = NCONS, nGradEqn = NGRAD, &
                                                  f = f, &
                                                  EllipticFlux = ViscousFlux0D, &
                                                  QLeft = f % storage(1) % Q(:,i,j), &
                                                  QRight = f % storage(2) % Q(:,i,j), &
                                                  U_xLeft = f % storage(1) % U_x(:,i,j), &
                                                  U_yLeft = f % storage(1) % U_y(:,i,j), &
                                                  U_zLeft = f % storage(1) % U_z(:,i,j), &
                                                  U_xRight = f % storage(2) % U_x(:,i,j), &
                                                  U_yRight = f % storage(2) % U_y(:,i,j), &
                                                  U_zRight = f % storage(2) % U_z(:,i,j), &
                                                  nHat = f % geom % normal(:,i,j) , &
                                                  dWall = f % geom % dWall(i,j), &
                                                  flux  = visc_flux(:,i,j) )

            end do
         end do

         else
         DO j = 0, f % Nf(2)
            DO i = 0, f % Nf(1)

!      
!              --------------
!              Viscous fluxes
!              --------------
!      
               CALL EllipticDiscretization % RiemannSolverWithSGS(f = f, &
                                                  QLeft = f % storage(1) % Q(:,i,j), &
                                                  QRight = f % storage(2) % Q(:,i,j), &
                                                  U_xLeft = f % storage(1) % U_x(:,i,j), &
                                                  U_yLeft = f % storage(1) % U_y(:,i,j), &
                                                  U_zLeft = f % storage(1) % U_z(:,i,j), &
                                                  U_xRight = f % storage(2) % U_x(:,i,j), &
                                                  U_yRight = f % storage(2) % U_y(:,i,j), &
                                                  U_zRight = f % storage(2) % U_z(:,i,j), &
                                                  nHat = f % geom % normal(:,i,j) , &
                                                  dWall = f % geom % dWall(i,j), &
                                                  flux  = visc_flux(:,i,j) )

            end do
         end do
         end if

         DO j = 0, f % Nf(2)
            DO i = 0, f % Nf(1)
!      
!              --------------
!              Invscid fluxes
!              --------------
!      
               CALL RiemannSolver(QLeft  = f % storage(1) % Q(:,i,j), &
                                  QRight = f % storage(2) % Q(:,i,j), &
                                  nHat   = f % geom % normal(:,i,j), &
                                  t1     = f % geom % t1(:,i,j), &
                                  t2     = f % geom % t2(:,i,j), &
                                  flux   = inv_flux(:,i,j) )

               
!
!              Multiply by the Jacobian
!              ------------------------
               flux(:,i,j) = ( inv_flux(:,i,j) - visc_flux(:,i,j)) * f % geom % jacobian(i,j)
               
            END DO   
         END DO  
!
!        ---------------------------
!        Return the flux to elements
!        ---------------------------
!
         call f % ProjectFluxToElements(NCONS, flux, (/1,2/))

      END SUBROUTINE computeElementInterfaceFlux_NS

      SUBROUTINE computeMPIFaceFlux_NS(f)
         use FaceClass
         use RiemannSolvers_NS
         IMPLICIT NONE
         TYPE(Face)   , INTENT(inout) :: f   
         integer       :: i, j
         integer       :: thisSide
         real(kind=RP) :: inv_flux(1:NCONS,0:f % Nf(1),0:f % Nf(2))
         real(kind=RP) :: visc_flux(1:NCONS,0:f % Nf(1),0:f % Nf(2))
         real(kind=RP) :: flux(1:NCONS,0:f % Nf(1),0:f % Nf(2))
!
!        --------------
!        Invscid fluxes
!        --------------
!
         DO j = 0, f % Nf(2)
            DO i = 0, f % Nf(1)
!      
!              --------------
!              Viscous fluxes
!              --------------
!      
               CALL EllipticDiscretization % RiemannSolver(nEqn = NCONS, nGradEqn = NGRAD, &
                                                  f = f, &
                                                  EllipticFlux = ViscousFlux0D, &
                                                  QLeft = f % storage(1) % Q(:,i,j), &
                                                  QRight = f % storage(2) % Q(:,i,j), &
                                                  U_xLeft = f % storage(1) % U_x(:,i,j), &
                                                  U_yLeft = f % storage(1) % U_y(:,i,j), &
                                                  U_zLeft = f % storage(1) % U_z(:,i,j), &
                                                  U_xRight = f % storage(2) % U_x(:,i,j), &
                                                  U_yRight = f % storage(2) % U_y(:,i,j), &
                                                  U_zRight = f % storage(2) % U_z(:,i,j), &
                                                  nHat = f % geom % normal(:,i,j) , &
                                                  dWall = f % geom % dWall(i,j), &
                                                  flux  = visc_flux(:,i,j) )

!
               CALL RiemannSolver(QLeft  = f % storage(1) % Q(:,i,j), &
                                  QRight = f % storage(2) % Q(:,i,j), &
                                  nHat   = f % geom % normal(:,i,j), &
                                  t1     = f % geom % t1(:,i,j), &
                                  t2     = f % geom % t2(:,i,j), &
                                  flux   = inv_flux(:,i,j) )
!
!              Multiply by the Jacobian
!              ------------------------
               flux(:,i,j) = ( inv_flux(:,i,j) - visc_flux(:,i,j)) * f % geom % jacobian(i,j)
               
            END DO   
         END DO  
!
!        ---------------------------
!        Return the flux to elements: The sign in eR % storage % FstarB has already been accouted.
!        ---------------------------
!
         thisSide = maxloc(f % elementIDs, dim = 1)
         call f % ProjectFluxToElements(NCONS, flux, (/thisSide, HMESH_NONE/))

      end subroutine ComputeMPIFaceFlux_NS

      SUBROUTINE computeBoundaryFlux_NS(f, time, externalStateProcedure , externalGradientsProcedure)
      USE ElementClass
      use FaceClass
      USE EllipticDiscretizations
      USE RiemannSolvers_NS
      USE BoundaryConditionFunctions
      IMPLICIT NONE
!
!     ---------
!     Arguments
!     ---------
!
      type(Face),    intent(inout) :: f
      REAL(KIND=RP)                :: time
      procedure(BCState_FCN)       :: externalState
      procedure(BCGradients_FCN)   :: externalGradients
!
!     ---------------
!     Local variables
!     ---------------
!
      INTEGER                         :: i, j
      INTEGER, DIMENSION(2)           :: N
      REAL(KIND=RP)                   :: inv_flux(NCONS)
      REAL(KIND=RP)                   :: UGradExt(NDIM , NGRAD) 
      real(kind=RP)                   :: visc_flux(NCONS, 0:f % Nf(1), 0:f % Nf(2))
      real(kind=RP)                   :: fStar(NCONS, 0:f % Nf(1), 0: f % Nf(2))
      
      CHARACTER(LEN=BC_STRING_LENGTH) :: boundaryType
      CHARACTER(LEN=BC_STRING_LENGTH) :: boundaryName
            
      boundaryType = f % boundaryType
      boundaryName = f % boundaryName
!
!     -------------------
!     Get external states
!     -------------------
!
      do j = 0, f % Nf(2)  ; do i = 0, f % Nf(1)
         f % storage(2) % Q(:,i,j) = f % storage(1) % Q(:,i,j)
         CALL externalStateProcedure( f % geom % x(:,i,j), &
                                      time, &
                                      f % geom % normal(:,i,j), &
                                      f % storage(2) % Q(:,i,j),&
                                      boundaryType, boundaryName )

      end do               ; end do

      if ( flowIsNavierStokes ) then
         do j = 0, f % Nf(2)  ; do i = 0, f % Nf(1)
            UGradExt(IX,:) = f % storage(1) % U_x(:,i,j)
            UGradExt(IY,:) = f % storage(1) % U_y(:,i,j)
            UGradExt(IZ,:) = f % storage(1) % U_z(:,i,j)

            CALL externalGradientsProcedure(  f % geom % x(:,i,j), &
                                              time, &
                                              f % geom % normal(:,i,j), &
                                              UGradExt,&
                                              boundaryType, boundaryName)

            f % storage(2) % U_x(:,i,j) = UGradExt(IX,:)
            f % storage(2) % U_y(:,i,j) = UGradExt(IY,:)
            f % storage(2) % U_z(:,i,j) = UGradExt(IZ,:)

!   
!           --------------
!           Viscous fluxes
!           --------------
!   
            CALL EllipticDiscretization % RiemannSolver(nEqn = NCONS, nGradEqn = NGRAD, &
                                               f = f, &
                                               EllipticFlux = ViscousFlux0D, &
                                               QLeft = f % storage(1) % Q(:,i,j), &
                                               QRight = f % storage(2) % Q(:,i,j), &
                                               U_xLeft = f % storage(1) % U_x(:,i,j), &
                                               U_yLeft = f % storage(1) % U_y(:,i,j), &
                                               U_zLeft = f % storage(1) % U_z(:,i,j), &
                                               U_xRight = f % storage(2) % U_x(:,i,j), &
                                               U_yRight = f % storage(2) % U_y(:,i,j), &
                                               U_zRight = f % storage(2) % U_z(:,i,j), &
                                               nHat = f % geom % normal(:,i,j) , &
                                               dWall = f % geom % dWall(i,j), &
                                               flux  = visc_flux(:,i,j) )

         end do               ; end do
      else
         visc_flux = 0.0_RP

      end if

      DO j = 0, f % Nf(2)
         DO i = 0, f % Nf(1)
!
!           Hyperbolic part
!           -------------
            CALL RiemannSolver(QLeft  = f % storage(1) % Q(:,i,j), &
                               QRight = f % storage(2) % Q(:,i,j), &   
                               nHat   = f % geom % normal(:,i,j), &
                               t1     = f % geom % t1(:,i,j), &
                               t2     = f % geom % t2(:,i,j), &
                               flux   = inv_flux)

            fStar(:,i,j) = (inv_flux - visc_flux(:,i,j) ) * f % geom % jacobian(i,j)
         END DO   
      END DO   

      call f % ProjectFluxToElements(NCONS, fStar, (/1, HMESH_NONE/))

      END SUBROUTINE computeBoundaryFlux_NS

      SUBROUTINE computeElementInterfaceFlux_SVV(f)
         use FaceClass
         use RiemannSolvers_NS
         IMPLICIT NONE
         TYPE(Face)   , INTENT(inout) :: f   
         integer       :: i, j
         real(kind=RP) :: inv_flux(1:NCONS,0:f % Nf(1),0:f % Nf(2))
         real(kind=RP) :: visc_flux(1:NCONS,0:f % Nf(1),0:f % Nf(2))
         real(kind=RP) :: SVV_flux(1:NCONS,0:f % Nf(1),0:f % Nf(2))
         real(kind=RP) :: flux(1:NCONS,0:f % Nf(1),0:f % Nf(2))
!
!        ----------
!        SVV fluxes
!        ----------
!
         if ( SVV % enabled ) then 
         CALL SVV % RiemannSolver(f = f, &
                        EllipticFlux = ViscousFlux2D, &
                              QLeft = f % storage(1) % Q, &
                             QRight = f % storage(2) % Q, &
                            U_xLeft = f % storage(1) % U_x, &
                            U_yLeft = f % storage(1) % U_y, &
                            U_zLeft = f % storage(1) % U_z, &
                           U_xRight = f % storage(2) % U_x, &
                           U_yRight = f % storage(2) % U_y, &
                           U_zRight = f % storage(2) % U_z, &
                              flux  = SVV_flux               )
         else
            SVV_Flux = 0.0_RP

         end if
!
         DO j = 0, f % Nf(2)
            DO i = 0, f % Nf(1)

!      
!              --------------
!              Viscous fluxes
!              --------------
!      
               CALL EllipticDiscretization % RiemannSolver(nEqn = NCONS, nGradEqn = NGRAD, &
                                                  f = f, &
                                                  EllipticFlux = ViscousFlux0D, &
                                                  QLeft = f % storage(1) % Q(:,i,j), &
                                                  QRight = f % storage(2) % Q(:,i,j), &
                                                  U_xLeft = f % storage(1) % U_x(:,i,j), &
                                                  U_yLeft = f % storage(1) % U_y(:,i,j), &
                                                  U_zLeft = f % storage(1) % U_z(:,i,j), &
                                                  U_xRight = f % storage(2) % U_x(:,i,j), &
                                                  U_yRight = f % storage(2) % U_y(:,i,j), &
                                                  U_zRight = f % storage(2) % U_z(:,i,j), &
                                                  nHat = f % geom % normal(:,i,j) , &
                                                  dWall = f % geom % dWall(i,j), &
                                                  flux  = visc_flux(:,i,j) )

!      
!              --------------
!              Invscid fluxes
!              --------------
!      
               CALL RiemannSolver(QLeft  = f % storage(1) % Q(:,i,j), &
                                  QRight = f % storage(2) % Q(:,i,j), &
                                  nHat   = f % geom % normal(:,i,j), &
                                  t1     = f % geom % t1(:,i,j), &
                                  t2     = f % geom % t2(:,i,j), &
                                  flux   = inv_flux(:,i,j) )

               
!
!              Multiply by the Jacobian
!              ------------------------
               flux(:,i,j) = ( inv_flux(:,i,j) - visc_flux(:,i,j) - SVV_flux(:,i,j) ) * f % geom % jacobian(i,j)
               
            END DO   
         END DO  
!
!        ---------------------------
!        Return the flux to elements
!        ---------------------------
!
         call f % ProjectFluxToElements(NCONS, flux, (/1,2/))

      END SUBROUTINE computeElementInterfaceFlux_SVV

      SUBROUTINE computeMPIFaceFlux_SVV(f)
         use FaceClass
         use RiemannSolvers_NS
         IMPLICIT NONE
         TYPE(Face)   , INTENT(inout) :: f   
         integer       :: i, j
         integer       :: thisSide
         real(kind=RP) :: inv_flux(1:NCONS,0:f % Nf(1),0:f % Nf(2))
         real(kind=RP) :: visc_flux(1:NCONS,0:f % Nf(1),0:f % Nf(2))
         real(kind=RP) :: SVV_flux(1:NCONS,0:f % Nf(1),0:f % Nf(2))
         real(kind=RP) :: flux(1:NCONS,0:f % Nf(1),0:f % Nf(2))
!
!        ----------
!        SVV fluxes
!        ----------
!
         CALL SVV % RiemannSolver(f = f, &
                            EllipticFlux = ViscousFlux2D, &
                              QLeft = f % storage(1) % Q, &
                             QRight = f % storage(2) % Q, &
                            U_xLeft = f % storage(1) % U_x, &
                            U_yLeft = f % storage(1) % U_y, &
                            U_zLeft = f % storage(1) % U_z, &
                           U_xRight = f % storage(2) % U_x, &
                           U_yRight = f % storage(2) % U_y, &
                           U_zRight = f % storage(2) % U_z, &
                              flux  = SVV_flux               )

!
!        --------------
!        Invscid fluxes
!        --------------
!
         DO j = 0, f % Nf(2)
            DO i = 0, f % Nf(1)
!      
!              --------------
!              Viscous fluxes
!              --------------
!      
               CALL EllipticDiscretization % RiemannSolver(nEqn = NCONS, nGradEqn = NGRAD, &
                                                  f = f, &
                                                  EllipticFlux = ViscousFlux0D, &
                                                  QLeft = f % storage(1) % Q(:,i,j), &
                                                  QRight = f % storage(2) % Q(:,i,j), &
                                                  U_xLeft = f % storage(1) % U_x(:,i,j), &
                                                  U_yLeft = f % storage(1) % U_y(:,i,j), &
                                                  U_zLeft = f % storage(1) % U_z(:,i,j), &
                                                  U_xRight = f % storage(2) % U_x(:,i,j), &
                                                  U_yRight = f % storage(2) % U_y(:,i,j), &
                                                  U_zRight = f % storage(2) % U_z(:,i,j), &
                                                  nHat = f % geom % normal(:,i,j) , &
                                                  dWall = f % geom % dWall(i,j), &
                                                  flux  = visc_flux(:,i,j) )

!
               CALL RiemannSolver(QLeft  = f % storage(1) % Q(:,i,j), &
                                  QRight = f % storage(2) % Q(:,i,j), &
                                  nHat   = f % geom % normal(:,i,j), &
                                  t1     = f % geom % t1(:,i,j), &
                                  t2     = f % geom % t2(:,i,j), &
                                  flux   = inv_flux(:,i,j) )
!
!              Multiply by the Jacobian
!              ------------------------
               flux(:,i,j) = ( inv_flux(:,i,j) - visc_flux(:,i,j) - SVV_flux(:,i,j) ) * f % geom % jacobian(i,j)
               
            END DO   
         END DO  
!
!        ---------------------------
!        Return the flux to elements: The sign in eR % storage % FstarB has already been accouted.
!        ---------------------------
!
         thisSide = maxloc(f % elementIDs, dim = 1)
         call f % ProjectFluxToElements(NCONS, flux, (/thisSide, HMESH_NONE/))


      end subroutine ComputeMPIFaceFlux_SVV

      SUBROUTINE computeBoundaryFlux_SVV(f, time, externalStateProcedure , externalGradientsProcedure)
      USE ElementClass
      use FaceClass
      USE EllipticDiscretizations
      USE RiemannSolvers_NS
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
      REAL(KIND=RP)                   :: inv_flux(NCONS)
      REAL(KIND=RP)                   :: UGradExt(NDIM , NGRAD) 
      real(kind=RP)                   :: visc_flux(NCONS, 0:f % Nf(1), 0:f % Nf(2))
      real(kind=RP)                   :: SVV_flux(NCONS, 0:f % Nf(1), 0:f % Nf(2))
      real(kind=RP)                   :: fStar(NCONS, 0:f % Nf(1), 0: f % Nf(2))
      
      CHARACTER(LEN=BC_STRING_LENGTH) :: boundaryType
            
      boundaryType = f % boundaryType
!
!     -------------------
!     Get external states
!     -------------------
!
      do j = 0, f % Nf(2)  ; do i = 0, f % Nf(1)
         f % storage(2) % Q(:,i,j) = f % storage(1) % Q(:,i,j)
         CALL externalStateProcedure( f % geom % x(:,i,j), &
                                      time, &
                                      f % geom % normal(:,i,j), &
                                      f % storage(2) % Q(:,i,j),&
                                      boundaryType )

      end do               ; end do

      if ( flowIsNavierStokes ) then
         do j = 0, f % Nf(2)  ; do i = 0, f % Nf(1)
            UGradExt(IX,:) = f % storage(1) % U_x(:,i,j)
            UGradExt(IY,:) = f % storage(1) % U_y(:,i,j)
            UGradExt(IZ,:) = f % storage(1) % U_z(:,i,j)

            CALL externalGradientsProcedure(  f % geom % x(:,i,j), &
                                              time, &
                                              f % geom % normal(:,i,j), &
                                              UGradExt,&
                                              boundaryType )

            f % storage(2) % U_x(:,i,j) = UGradExt(IX,:)
            f % storage(2) % U_y(:,i,j) = UGradExt(IY,:)
            f % storage(2) % U_z(:,i,j) = UGradExt(IZ,:)

!   
!           --------------
!           Viscous fluxes
!           --------------
!   
            CALL EllipticDiscretization % RiemannSolver(nEqn = NCONS, nGradEqn = NGRAD, &
                                               f = f, &
                                               EllipticFlux = ViscousFlux0D, &
                                               QLeft = f % storage(1) % Q(:,i,j), &
                                               QRight = f % storage(2) % Q(:,i,j), &
                                               U_xLeft = f % storage(1) % U_x(:,i,j), &
                                               U_yLeft = f % storage(1) % U_y(:,i,j), &
                                               U_zLeft = f % storage(1) % U_z(:,i,j), &
                                               U_xRight = f % storage(2) % U_x(:,i,j), &
                                               U_yRight = f % storage(2) % U_y(:,i,j), &
                                               U_zRight = f % storage(2) % U_z(:,i,j), &
                                               nHat = f % geom % normal(:,i,j) , &
                                               dWall = f % geom % dWall(i,j), &
                                               flux  = visc_flux(:,i,j) )

         end do               ; end do
      else
         visc_flux = 0.0_RP

      end if
!
!     ----------
!     SVV fluxes
!     ----------
!
      CALL SVV % RiemannSolver(f = f, &
                           EllipticFlux = ViscousFlux2D, &
                           QLeft = f % storage(1) % Q, &
                          QRight = f % storage(2) % Q, &
                         U_xLeft = f % storage(1) % U_x, &
                         U_yLeft = f % storage(1) % U_y, &
                         U_zLeft = f % storage(1) % U_z, &
                        U_xRight = f % storage(2) % U_x, &
                        U_yRight = f % storage(2) % U_y, &
                        U_zRight = f % storage(2) % U_z, &
                           flux  = SVV_flux               )

      DO j = 0, f % Nf(2)
         DO i = 0, f % Nf(1)
!
!           Hyperbolic part
!           -------------
            CALL RiemannSolver(QLeft  = f % storage(1) % Q(:,i,j), &
                               QRight = f % storage(2) % Q(:,i,j), &   
                               nHat   = f % geom % normal(:,i,j), &
                               t1     = f % geom % t1(:,i,j), &
                               t2     = f % geom % t2(:,i,j), &
                               flux   = inv_flux)

            fStar(:,i,j) = (inv_flux - visc_flux(:,i,j) - SVV_flux(:,i,j) ) * f % geom % jacobian(i,j)
         END DO   
      END DO   

      call f % ProjectFluxToElements(NCONS, fStar, (/1, HMESH_NONE/))

      END SUBROUTINE computeBoundaryFlux_SVV
!
!////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
!
!     -----------------------------------------------------------------
!     Subroutine to compute the isolated fluxes on both sides of a face
!     -----------------------------------------------------------------
      SUBROUTINE computeIsolatedFaceFluxes_NS(f)
         use FaceClass
         IMPLICIT NONE
         TYPE(Face)   , INTENT(inout) :: f   
         !-----------------------------------------------------------
         integer       :: i, j
         real(kind=RP) :: inv_fluxL (1:NCONS,0:f % Nf(1),0:f % Nf(2))
         real(kind=RP) :: inv_fluxR (1:NCONS,0:f % Nf(1),0:f % Nf(2))
         real(kind=RP) :: visc_fluxL(1:NCONS,0:f % Nf(1),0:f % Nf(2))
         real(kind=RP) :: visc_fluxR(1:NCONS,0:f % Nf(1),0:f % Nf(2))
         real(kind=RP) :: fluxL     (1:NCONS,0:f % Nf(1),0:f % Nf(2))
         real(kind=RP) :: fluxR     (1:NCONS,0:f % Nf(1),0:f % Nf(2))
         !-----------------------------------------------------------
         real(kind=RP) :: nHat(NDIM)
         real(kind=RP) :: mu, kappa             !
         real(kind=RP) :: flux_vec(NCONS,NDIM)  ! Flux tensor
         !-----------------------------------------------------------
         
#if defined(NAVIERSTOKES)
         mu    = dimensionless % mu
         kappa = dimensionless % kappa
#else
         mu = 0.0_RP
         kappa = 0.0_RP
#endif
         
         if ( .not. LESModel % active ) then
            do j = 0, f % Nf(2)
               do i = 0, f % Nf(1)
                  
                  nHat = f % geom % normal (:,i,j)
                  
!                 
!                 Viscous flux on the left
!                 ------------------------
                  
                  call ViscousFlux(nEqn  = NCONS,                       &
                                   nGradEqn = NGRAD,                    &
                                   Q     = f % storage(1) % Q  (:,i,j), &
                                   U_x   = f % storage(1) % U_x(:,i,j), &
                                   U_y   = f % storage(1) % U_y(:,i,j), &
                                   U_z   = f % storage(1) % U_z(:,i,j), &
                                   mu    = mu                         , &
                                   kappa = kappa                      , &
                                   F     = flux_vec                     )
                  
                  visc_fluxL(:,i,j) = flux_vec(:,IX) * nHat(IX) + flux_vec(:,IY) * nHat(IY) + flux_vec(:,IZ) * nHat(IZ)
                  
!                 
!                 Viscous flux on the right
!                 -------------------------
                  
                  call ViscousFlux(nEqn  = NCONS, nGradEqn = NGRAD,     &
                                   Q     = f % storage(2) % Q  (:,i,j), &
                                   U_x   = f % storage(2) % U_x(:,i,j), &
                                   U_y   = f % storage(2) % U_y(:,i,j), &
                                   U_z   = f % storage(2) % U_z(:,i,j), &
                                   mu    = mu                         , &
                                   kappa = kappa                      , &
                                   F     = flux_vec                     )
                  
                  visc_fluxR(:,i,j) = flux_vec(:,IX) * nHat(IX) + flux_vec(:,IY) * nHat(IY) + flux_vec(:,IZ) * nHat(IZ) ! The inversion is performed later
               end do
            end do

         else
            error stop ':: The isolated face flux computation is not yet implemented for LES models. Sorry...'
         end if

         do j = 0, f % Nf(2)
            do i = 0, f % Nf(1)
               
               nHat = f % geom % normal (:,i,j)
               
!      
!              Hyperbolic flux on the left
!              -------------------------
               
               call EulerFlux(f % storage(1) % Q  (:,i,j), flux_vec)
               
               inv_fluxL(:,i,j) = flux_vec(:,IX) * nHat(IX) + flux_vec(:,IY) * nHat(IY) + flux_vec(:,IZ) * nHat(IZ)
               
!      
!              Hyperbolic flux on the right
!              --------------------------
               
               call EulerFlux(f % storage(2) % Q  (:,i,j), flux_vec)
               
               inv_fluxR(:,i,j) = flux_vec(:,IX) * nHat(IX) + flux_vec(:,IY) * nHat(IY) + flux_vec(:,IZ) * nHat(IZ) ! The inversion is performed later
               
!
!              Multiply by the Jacobian
!              ------------------------
               fluxL(:,i,j) = ( inv_fluxL(:,i,j) - visc_fluxL(:,i,j)) * f % geom % jacobian(i,j)
               
               fluxR(:,i,j) = ( inv_fluxR(:,i,j) - visc_fluxR(:,i,j)) * f % geom % jacobian(i,j)
               
            end do
         end do
!
!        ---------------------------
!        Return the flux to elements
!        ---------------------------
!
         ! Left element
         call f % ProjectFluxToElements(NCONS, fluxL, (/1,HMESH_NONE/))
         
         ! Right element
         call f % ProjectFluxToElements(NCONS, fluxR, (/2,HMESH_NONE/))

      END SUBROUTINE computeIsolatedFaceFluxes_NS
!
!////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
!
!     --------------------------------------------------------------------
!     Subroutine to compute the isolated fluxes on only one side of a face
!     --------------------------------------------------------------------
      SUBROUTINE computeIsolatedFaceFlux_NS(f,thisSide)
         use FaceClass
         IMPLICIT NONE
         !-----------------------------------------------------------
         TYPE(Face)   , INTENT(inout) :: f
         integer      , intent(in)    :: thisSide
         !-----------------------------------------------------------
         integer       :: i, j
         real(kind=RP) :: inv_flux (1:NCONS,0:f % Nf(1),0:f % Nf(2))
         real(kind=RP) :: visc_flux(1:NCONS,0:f % Nf(1),0:f % Nf(2))
         real(kind=RP) :: flux     (1:NCONS,0:f % Nf(1),0:f % Nf(2))    
         real(kind=RP) :: nHat(NDIM)                                    ! Face normal vector
         real(kind=RP) :: mu, kappa                                     ! dimensionless constants
         real(kind=RP) :: flux_vec(NCONS,NDIM)                          ! Flux tensor
         !-----------------------------------------------------------
         
#if defined(NAVIERSTOKES)
         mu    = dimensionless % mu
         kappa = dimensionless % kappa
#else
         mu = 0.0_RP
         kappa = 0.0_RP
#endif
         
         if ( .not. LESModel % active ) then
            do j = 0, f % Nf(2)
               do i = 0, f % Nf(1)
                  
                  nHat = f % geom % normal (:,i,j)
                  
!                 
!                 Viscous flux
!                 ------------
                  
                  call ViscousFlux(nEqn  = NCONS, nGradEqn = NGRAD,            &
                                   Q     = f % storage(thisSide) % Q  (:,i,j), &
                                   U_x   = f % storage(thisSide) % U_x(:,i,j), &
                                   U_y   = f % storage(thisSide) % U_y(:,i,j), &
                                   U_z   = f % storage(thisSide) % U_z(:,i,j), &
                                   mu    = mu                                , &
                                   kappa = kappa                             , &
                                   F     = flux_vec                          )
                  
                  visc_flux(:,i,j) = flux_vec(:,IX) * nHat(IX) + flux_vec(:,IY) * nHat(IY) + flux_vec(:,IZ) * nHat(IZ)
                  
               end do
            end do

         else
            error stop ':: The isolated face flux computation is not yet implemented for LES models. Sorry...'
         end if

         do j = 0, f % Nf(2)
            do i = 0, f % Nf(1)
               
               nHat = f % geom % normal (:,i,j)
               
!      
!              Hyperbolic flux
!              -------------
               
               call EulerFlux(f % storage(thisSide) % Q  (:,i,j), flux_vec)
               
               inv_flux(:,i,j) = flux_vec(:,IX) * nHat(IX) + flux_vec(:,IY) * nHat(IY) + flux_vec(:,IZ) * nHat(IZ) ! The inversion is performed later
               
!
!              Multiply by the Jacobian
!              ------------------------
               flux(:,i,j) = ( inv_flux(:,i,j) - visc_flux(:,i,j)) * f % geom % jacobian(i,j)
               
            end do
         end do
!
!        ------------------------------
!        Return the flux to the element
!        ------------------------------
!
         call f % ProjectFluxToElements(NCONS, flux, (/thisSide,HMESH_NONE/))

      END SUBROUTINE computeIsolatedFaceFlux_NS
!
!////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
!
!              GRADIENT PROCEDURES
!              -------------------
!
!////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
!
      subroutine DGSpatial_ComputeGradient( mesh , time , externalStateProcedure)
         use HexMeshClass
         implicit none
         type(HexMesh)                  :: mesh
         real(kind=RP),      intent(in) :: time
         procedure(BCState_FCN)         :: externalStateProcedure

         call EllipticDiscretization % ComputeGradient( NCONS, NGRAD, mesh , time , externalStateProcedure, NSGradientValuesForQ_0D, NSGradientValuesForQ_3D)

      end subroutine DGSpatial_ComputeGradient
!
!////////////////////////////////////////////////////////////////////////////////////////
!
end module SpatialDiscretization
