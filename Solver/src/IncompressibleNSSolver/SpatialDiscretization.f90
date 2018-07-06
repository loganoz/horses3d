!
!//////////////////////////////////////////////////////
!
!   @File:    SpatialDiscretization.f90
!   @Author:  Juan Manzanero (juan.manzanero@upm.es)
!   @Created: Wed Jun 20 18:14:45 2018
!   @Last revision date: Fri Jul  6 16:56:22 2018
!   @Last revision author: Juan Manzanero (juan.manzanero@upm.es)
!   @Last revision commit: 8261cb5f0eda85742f106865db8074403435942b
!
!//////////////////////////////////////////////////////
!
#include "Includes.h"
module SpatialDiscretization
      use SMConstants
      use HyperbolicDiscretizations
      use EllipticDiscretizations
      use DGIntegrals
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
      use VariableConversion, only: iNSGradientValuesForQ_0D, iNSGradientValuesForQ_3D, GetiNSOneFluidViscosity, GetiNSTwoFluidsViscosity
      use ProblemFileFunctions
#ifdef _HAS_MPI_
      use mpi
#endif

      private
      public   ComputeTimeDerivative, ComputeTimeDerivativeIsolated, viscousDiscretizationKey
      public   Initialize_SpaceAndTimeMethods, Finalize_SpaceAndTimeMethods

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
            use FaceClass,  only: Face
            use DGSEMClass, only: BCState_FCN, BCGradients_FCN
            IMPLICIT NONE
            type(Face),    intent(inout) :: f
            REAL(KIND=RP)                :: time
            PROCEDURE(BCState_FCN)       :: externalStateProcedure
            PROCEDURE(BCGradients_FCN)   :: externalGradientsProcedure
         end subroutine computeBoundaryFluxF
      end interface
      
      character(len=LINE_LENGTH), parameter  :: viscousDiscretizationKey = "viscous discretization"
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
!
!        ---------------
!        Local variables
!        ---------------
!
         character(len=LINE_LENGTH)       :: inviscidDiscretizationName
         character(len=LINE_LENGTH)       :: viscousDiscretizationName
         
         if (.not. mesh % child) then ! If this is a child mesh, all these constructs were already initialized for the parent mesh
         
            if ( MPI_Process % isRoot ) then
               write(STD_OUT,'(/)')
               call Section_Header("Spatial discretization scheme")
               write(STD_OUT,'(/)')
            end if
   !
   !        Initialize inviscid discretization
   !        ----------------------------------
            inviscidDiscretizationName = controlVariables % stringValueForKey(inviscidDiscretizationKey,requestedLength = LINE_LENGTH)

            call toLower(inviscidDiscretizationName)
         
            select case ( trim(inviscidDiscretizationName) )
            case ( "standard" )
               if (.not. allocated(HyperbolicDiscretization)) allocate( StandardDG_t  :: HyperbolicDiscretization )

            case ( "split-form")
               if (.not. allocated(HyperbolicDiscretization)) allocate( SplitDG_t     :: HyperbolicDiscretization)

            case default
               write(STD_OUT,'(A,A,A)') 'Requested inviscid discretization "',trim(inviscidDiscretizationName),'" is not implemented.'
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
               if ( .not. controlVariables % ContainsKey(viscousDiscretizationKey) ) then
                  print*, "Input file is missing entry for keyword: viscous discretization"
                  errorMessage(STD_OUT)
                  stop
               end if

               viscousDiscretizationName = controlVariables % stringValueForKey(viscousDiscretizationKey, requestedLength = LINE_LENGTH)
               call toLower(viscousDiscretizationName)
               
               select case ( trim(viscousDiscretizationName) )
               case("br1")
                  allocate(BassiRebay1_t     :: ViscousDiscretization)

               case("br2")
                  allocate(BassiRebay2_t     :: ViscousDiscretization)

               case("ip")
                  allocate(InteriorPenalty_t :: ViscousDiscretization)

               case default
                  write(STD_OUT,'(A,A,A)') 'Requested viscous discretization "',trim(viscousDiscretizationName),'" is not implemented.'
                  write(STD_OUT,'(A)') "Implemented discretizations are:"
                  write(STD_OUT,'(A)') "  * BR1"
                  write(STD_OUT,'(A)') "  * BR2"
                  write(STD_OUT,'(A)') "  * IP"
                  errorMessage(STD_OUT)
                  stop 

               end select

               select case (thermodynamics % number_of_fluids)
               case(1)
                  call ViscousDiscretization % Construct(controlVariables, iViscousFlux0D, iViscousFlux2D, iViscousFlux3D, GetiNSOneFluidViscosity, "NS")
               case(2)
                  call ViscousDiscretization % Construct(controlVariables, iViscousFlux0D, iViscousFlux2D, iViscousFlux3D, GetiNSTwoFluidsViscosity, "NS")
               end select

               call ViscousDiscretization % Describe
!
!        Compute wall distances
!        ----------------------
         call mesh % ComputeWallDistances
         
         end if

      end subroutine Initialize_SpaceAndTimeMethods
!
!////////////////////////////////////////////////////////////////////////
!
      subroutine Finalize_SpaceAndTimeMethods
         implicit none
         IF ( ALLOCATED(HyperbolicDiscretization) ) DEALLOCATE( HyperbolicDiscretization )
      end subroutine Finalize_SpaceAndTimeMethods
!
!////////////////////////////////////////////////////////////////////////
!
      SUBROUTINE ComputeTimeDerivative( mesh, particles, time, BCFunctions, mode)
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
         integer,             intent(in) :: mode
!
!        ---------------
!        Local variables
!        ---------------
!
         INTEGER :: k
!
!        Apply a limiter to the density
!        ------------------------------
         if ( enableDensityLimiter ) then
            do k = 1, mesh % no_of_elements
               call DensityLimiter(mesh % elements(k) % Nxyz, mesh % elements(k) % storage % Q)
            end do
         end if
!
!        -----------------------------------------
!        Prolongation of the solution to the faces
!        -----------------------------------------
!
!$omp parallel shared(mesh, time)
         call mesh % ProlongSolutionToFaces(NINC)
!
!        ----------------
!        Update MPI Faces
!        ----------------
!
#ifdef _HAS_MPI_
!$omp single
         call mesh % UpdateMPIFacesSolution(NINC)
!$omp end single
#endif
!
!        -----------------
!        Compute gradients
!        -----------------
!
         if ( computeGradients ) then
            CALL DGSpatial_ComputeGradient(mesh , time , BCFunctions(1) % externalState)
         end if

#ifdef _HAS_MPI_
!$omp single
         call mesh % UpdateMPIFacesGradients(NINC)
!$omp end single
#endif
!
!        -----------------------
!        Compute time derivative
!        -----------------------
!
         call ComputeNSTimeDerivative(mesh = mesh , &
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
      SUBROUTINE ComputeTimeDerivativeIsolated( mesh, particles, time, BCFunctions, mode)
         use EllipticDiscretizationClass
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
         integer,             intent(in) :: mode
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
         call mesh % ProlongSolutionToFaces(NINC)
!
!        -----------------------------------------------------
!        Compute LOCAL gradients and prolong them to the faces
!        -----------------------------------------------------
!
         if ( computeGradients ) then
            CALL BaseClass_ComputeGradient( ViscousDiscretization, NINC, NINC, mesh , time , BCFunctions(1) % externalState, iNSGradientValuesForQ_0D, iNSGradientValuesForQ_3D )
!
!           The prolongation is usually done in the viscous methods, but not in the BaseClass
!           ---------------------------------------------------------------------------------
            call mesh % ProlongGradientsToFaces(NINC)
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
!
!////////////////////////////////////////////////////////////////////////////////////
!
!           Navier--Stokes procedures
!           -------------------------
!
!////////////////////////////////////////////////////////////////////////////////////
!
      subroutine ComputeNSTimeDerivative( mesh , particles, t, externalState, externalGradients )
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
               CALL computeElementInterfaceFlux_iNS( f ) 
 
            case (HMESH_BOUNDARY) 
               CALL computeBoundaryFlux_iNS(f, t, externalState, externalGradients) 
 
            end select 
            end associate 
         end do 
!$omp end do 
!
!        **************************************************************
!        Surface integrals and scaling of elements without shared faces
!        **************************************************************
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
            call mesh % GatherMPIFacesGradients(NINC)
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
                  CALL computeMPIFaceFlux_iNS( f ) 
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

         if (.not. mesh % child) then
            if ( particles % active ) then             
!$omp do schedule(runtime)
               do eID = 1, size(mesh % elements)
                  call particles % AddSource(mesh % elements(eID), t, thermodynamics, dimensionless, refValues)
               end do
!$omp end do
            endif 
         end if
!
!        ***********
!        Add gravity
!        ***********
!
!$omp do schedule(runtime) private(i,j,k)
            do eID = 1, size(mesh % elements)
               associate(e => mesh % elements(eID))
               do k = 0, e % Nxyz(3) ; do j = 0, e % Nxyz(2) ; do i = 0, e % Nxyz(1)
                  e % storage % QDot(INSU:INSW,i,j,k) = e % storage % QDot(INSU:INSW,i,j,k) + &
                                                        e % storage % Q(INSRHO,i,j,k) * &
                                    dimensionless % invFr2 * dimensionless % gravity_dir

               end do                ; end do                ; end do
               end associate
            end do
!$omp end do
!
!        ******************************************************
!        Apply the chain rule to get velocities time derivative
!        ******************************************************
!
!$omp do schedule(runtime) private(i,j,k)
         do eID = 1, size(mesh % elements) 
            associate(e => mesh % elements(eID)) 
            do k = 0, e % Nxyz(3) ; do j = 0, e % Nxyz(2) ; do i = 0, e % Nxyz(1) 
               e % storage % QDot(INSU,i,j,k) = (e % storage % QDot(INSU,i,j,k) - e % storage % Q(INSU,i,j,k)*e % storage % QDot(INSRHO,i,j,k)) / e % storage % Q(INSRHO,i,j,k)
               e % storage % QDot(INSV,i,j,k) = (e % storage % QDot(INSV,i,j,k) - e % storage % Q(INSV,i,j,k)*e % storage % QDot(INSRHO,i,j,k)) / e % storage % Q(INSRHO,i,j,k)
               e % storage % QDot(INSW,i,j,k) = (e % storage % QDot(INSW,i,j,k) - e % storage % Q(INSW,i,j,k)*e % storage % QDot(INSRHO,i,j,k)) / e % storage % Q(INSRHO,i,j,k)
            end do         ; end do          ; end do 
            end associate 
         end do
!$omp end do
!
!           ***************
!           Add source term
!           ***************
!$omp do schedule(runtime) private(i,j,k)
            do eID = 1, mesh % no_of_elements
               associate ( e => mesh % elements(eID) )
               do k = 0, e % Nxyz(3)   ; do j = 0, e % Nxyz(2) ; do i = 0, e % Nxyz(1)
                  call UserDefinedSourceTermNS(e % geom % x(:,i,j,k), e % storage % Q(:,i,j,k), t, e % storage % S_NS(:,i,j,k), thermodynamics, dimensionless, refValues)
               end do                  ; end do                ; end do
               end associate
            end do
!$omp end do

!$omp do schedule(runtime) private(i,j,k)
            do eID = 1, mesh % no_of_elements
               associate ( e => mesh % elements(eID) )
               do k = 0, e % Nxyz(3)   ; do j = 0, e % Nxyz(2) ; do i = 0, e % Nxyz(1)
                  e % storage % QDot(:,i,j,k) = e % storage % QDot(:,i,j,k) + e % storage % S_NS(:,i,j,k)
               end do                  ; end do                ; end do
               end associate
            end do
!$omp end do


      end subroutine ComputeNSTimeDerivative
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
!
!        ****************
!        Volume integrals
!        ****************
!
!$omp do schedule(runtime) 
         do eID = 1 , size(mesh % elements)
            call TimeDerivative_StrongVolumetricContribution( mesh % elements(eID) , t)
         end do
!$omp end do
!
!        *******************
!        Scaling of elements
!        *******************
! 
!$omp do schedule(runtime) private(i,j,k)
         do eID = 1, size(mesh % elements) 
            associate(e => mesh % elements(eID))

            do k = 0, e % Nxyz(3) ; do j = 0, e % Nxyz(2) ; do i = 0, e % Nxyz(1) 
               e % storage % QDot(:,i,j,k) = e % storage % QDot(:,i,j,k) / e % geom % jacobian(i,j,k) 
            end do         ; end do          ; end do 
            end associate 
         end do
!$omp end do
         
      end subroutine TimeDerivative_ComputeQDotIsolated
!
!///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
!
      subroutine TimeDerivative_StrongVolumetricContribution( e , t )
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
         real(kind=RP) :: inviscidContravariantFlux ( 1:NINC, 0:e%Nxyz(1) , 0:e%Nxyz(2) , 0:e%Nxyz(3), 1:NDIM ) 
         real(kind=RP) :: fSharp(1:NINC, 0:e%Nxyz(1), 0:e%Nxyz(1), 0:e%Nxyz(2), 0:e%Nxyz(3))
         real(kind=RP) :: gSharp(1:NINC, 0:e%Nxyz(2), 0:e%Nxyz(1), 0:e%Nxyz(2), 0:e%Nxyz(3))
         real(kind=RP) :: hSharp(1:NINC, 0:e%Nxyz(3), 0:e%Nxyz(1), 0:e%Nxyz(2), 0:e%Nxyz(3))
         real(kind=RP) :: viscousContravariantFlux  ( 1:NINC, 0:e%Nxyz(1) , 0:e%Nxyz(2) , 0:e%Nxyz(3), 1:NDIM ) 
         real(kind=RP) :: contravariantFlux         ( 1:NINC, 0:e%Nxyz(1) , 0:e%Nxyz(2) , 0:e%Nxyz(3), 1:NDIM ) 
         integer       :: eID
!
!        *************************************
!        Compute interior contravariant fluxes
!        *************************************
!
!        Compute inviscid contravariant flux
!        -----------------------------------
         call HyperbolicDiscretization % ComputeInnerFluxes ( e , iEulerFlux3D, inviscidContravariantFlux ) 
!
!        Compute viscous contravariant flux
!        ----------------------------------
         call ViscousDiscretization  % ComputeInnerFluxes ( NINC, NINC, e , viscousContravariantFlux) 
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
            contravariantFlux = inviscidContravariantFlux - viscousContravariantFlux 
!
!           Perform the Weak Volume Green integral
!           --------------------------------------
            e % storage % QDot = ScalarStrongIntegrals % StdVolumeGreen ( e , NINC, contravariantFlux ) 

         type is (SplitDG_t)
            ERROR stop ':: TimeDerivative_StrongVolumetricContribution not implemented for split form'
!~ !
!~ !           Compute sharp fluxes for skew-symmetric approximations
!~ !           ------------------------------------------------------
!~             call HyperbolicDiscretization % ComputeSplitFormFluxes(e, inviscidContravariantFlux, fSharp, gSharp, hSharp)
!~ !
!~ !           Peform the Weak volume green integral
!~ !           -------------------------------------
!~             viscousContravariantFlux = viscousContravariantFlux + SVVContravariantFlux

!~             e % storage % QDot = -ScalarWeakIntegrals % SplitVolumeDivergence( e, fSharp, gSharp, hSharp, viscousContravariantFlux)

         end select

      end subroutine TimeDerivative_StrongVolumetricContribution
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
         real(kind=RP) :: inviscidContravariantFlux ( 1:NINC, 0:e%Nxyz(1) , 0:e%Nxyz(2) , 0:e%Nxyz(3), 1:NDIM ) 
         real(kind=RP) :: fSharp(1:NINC, 0:e%Nxyz(1), 0:e%Nxyz(1), 0:e%Nxyz(2), 0:e%Nxyz(3))
         real(kind=RP) :: gSharp(1:NINC, 0:e%Nxyz(2), 0:e%Nxyz(1), 0:e%Nxyz(2), 0:e%Nxyz(3))
         real(kind=RP) :: hSharp(1:NINC, 0:e%Nxyz(3), 0:e%Nxyz(1), 0:e%Nxyz(2), 0:e%Nxyz(3))
         real(kind=RP) :: viscousContravariantFlux  ( 1:NINC, 0:e%Nxyz(1) , 0:e%Nxyz(2) , 0:e%Nxyz(3), 1:NDIM ) 
         real(kind=RP) :: contravariantFlux         ( 1:NINC, 0:e%Nxyz(1) , 0:e%Nxyz(2) , 0:e%Nxyz(3), 1:NDIM ) 
         integer       :: eID
!
!        *************************************
!        Compute interior contravariant fluxes
!        *************************************
!
!        Compute inviscid contravariant flux
!        -----------------------------------
         call HyperbolicDiscretization % ComputeInnerFluxes ( e , iEulerFlux3D, inviscidContravariantFlux ) 
!
!        Compute viscous contravariant flux
!        ----------------------------------
         call ViscousDiscretization  % ComputeInnerFluxes ( NINC, NINC, e , viscousContravariantFlux) 
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
            contravariantFlux = inviscidContravariantFlux - viscousContravariantFlux 
!
!           Perform the Weak Volume Green integral
!           --------------------------------------
            e % storage % QDot = ScalarWeakIntegrals % StdVolumeGreen ( e, NINC, contravariantFlux ) 

         type is (SplitDG_t)
!
!           Compute sharp fluxes for skew-symmetric approximations
!           ------------------------------------------------------
            call HyperbolicDiscretization % ComputeSplitFormFluxes(e, inviscidContravariantFlux, fSharp, gSharp, hSharp)
!
!           Peform the Weak volume green integral
!           -------------------------------------
            viscousContravariantFlux = viscousContravariantFlux 

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

         e % storage % QDot = e % storage % QDot - ScalarWeakIntegrals % StdFace( e, NINC, &
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
      SUBROUTINE computeElementInterfaceFlux_iNS(f)
         use FaceClass
         use RiemannSolvers_iNS
         IMPLICIT NONE
         TYPE(Face)   , INTENT(inout) :: f   
         integer       :: i, j
         real(kind=RP) :: inv_flux(1:NINC,0:f % Nf(1),0:f % Nf(2))
         real(kind=RP) :: visc_flux(1:NINC,0:f % Nf(1),0:f % Nf(2))
         real(kind=RP) :: flux(1:NINC,0:f % Nf(1),0:f % Nf(2))
         real(kind=RP) :: muL, muR, mu

         DO j = 0, f % Nf(2)
            DO i = 0, f % Nf(1)

               call ViscousDiscretization % GetViscosity(f % storage(1) % Q(INSRHO,i,j), muL)
               call ViscousDiscretization % GetViscosity(f % storage(2) % Q(INSRHO,i,j), muR)
               mu = 0.5_RP * (muL + muR)
!      
!              --------------
!              Viscous fluxes
!              --------------
!      
               CALL ViscousDiscretization % RiemannSolver(nEqn = NINC, nGradEqn = NINC, &
                                                  f = f, &
                                                  QLeft = f % storage(1) % Q(:,i,j), &
                                                  QRight = f % storage(2) % Q(:,i,j), &
                                                  U_xLeft = f % storage(1) % U_x(:,i,j), &
                                                  U_yLeft = f % storage(1) % U_y(:,i,j), &
                                                  U_zLeft = f % storage(1) % U_z(:,i,j), &
                                                  U_xRight = f % storage(2) % U_x(:,i,j), &
                                                  U_yRight = f % storage(2) % U_y(:,i,j), &
                                                  U_zRight = f % storage(2) % U_z(:,i,j), &
                                                  mu   = mu, &
                                                  nHat = f % geom % normal(:,i,j) , &
                                                  dWall = f % geom % dWall(i,j), &
                                                  flux  = visc_flux(:,i,j) )

            end do
         end do

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
         call f % ProjectFluxToElements(NINC, flux, (/1,2/))

      END SUBROUTINE computeElementInterfaceFlux_iNS

      SUBROUTINE computeMPIFaceFlux_iNS(f)
         use FaceClass
         use RiemannSolvers_iNS
         IMPLICIT NONE
         TYPE(Face)   , INTENT(inout) :: f   
         integer       :: i, j
         integer       :: thisSide
         real(kind=RP) :: inv_flux(1:NINC,0:f % Nf(1),0:f % Nf(2))
         real(kind=RP) :: visc_flux(1:NINC,0:f % Nf(1),0:f % Nf(2))
         real(kind=RP) :: flux(1:NINC,0:f % Nf(1),0:f % Nf(2))
         real(kind=RP) :: mu
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
               call ViscousDiscretization % GetViscosity(f % storage(1) % Q(INSRHO,i,j), mu)

               CALL ViscousDiscretization % RiemannSolver(nEqn = NINC, nGradEqn = NINC, &
                                                  f = f, &
                                                  QLeft = f % storage(1) % Q(:,i,j), &
                                                  QRight = f % storage(2) % Q(:,i,j), &
                                                  U_xLeft = f % storage(1) % U_x(:,i,j), &
                                                  U_yLeft = f % storage(1) % U_y(:,i,j), &
                                                  U_zLeft = f % storage(1) % U_z(:,i,j), &
                                                  U_xRight = f % storage(2) % U_x(:,i,j), &
                                                  U_yRight = f % storage(2) % U_y(:,i,j), &
                                                  U_zRight = f % storage(2) % U_z(:,i,j), &
                                                  mu   = mu, &
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
         call f % ProjectFluxToElements(NINC, flux, (/thisSide, HMESH_NONE/))

      end subroutine ComputeMPIFaceFlux_iNS

      SUBROUTINE computeBoundaryFlux_iNS(f, time, externalStateProcedure , externalGradientsProcedure)
      USE ElementClass
      use FaceClass
      USE RiemannSolvers_iNS
      USE BoundaryConditionFunctions
      IMPLICIT NONE
!
!     ---------
!     Arguments
!     ---------
!
      type(Face),    intent(inout) :: f
      REAL(KIND=RP)                :: time
      procedure(BCState_FCN)       :: externalStateProcedure
      procedure(BCGradients_FCN)   :: externalGradientsProcedure
!
!     ---------------
!     Local variables
!     ---------------
!
      INTEGER                         :: i, j
      INTEGER, DIMENSION(2)           :: N
      REAL(KIND=RP)                   :: inv_flux(NINC)
      REAL(KIND=RP)                   :: UGradExt(NDIM , NINC) 
      real(kind=RP)                   :: visc_flux(NINC, 0:f % Nf(1), 0:f % Nf(2))
      real(kind=RP)                   :: fStar(NINC, 0:f % Nf(1), 0: f % Nf(2))
      real(kind=RP)                   :: mu
      
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
         CALL externalStateProcedure( NINC, &
                                      f % geom % x(:,i,j), &
                                      time, &
                                      f % geom % normal(:,i,j), &
                                      f % storage(2) % Q(:,i,j),&
                                      boundaryType, boundaryName )

      end do               ; end do

      if ( .true. ) then
         do j = 0, f % Nf(2)  ; do i = 0, f % Nf(1)
            UGradExt(IX,:) = f % storage(1) % U_x(:,i,j)
            UGradExt(IY,:) = f % storage(1) % U_y(:,i,j)
            UGradExt(IZ,:) = f % storage(1) % U_z(:,i,j)

            CALL externalGradientsProcedure(  NINC, &
                                              f % geom % x(:,i,j), &
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
            call ViscousDiscretization % GetViscosity(f % storage(1) % Q(INSRHO,i,j), mu)

            CALL ViscousDiscretization % RiemannSolver(nEqn = NINC, nGradEqn = NINC, &
                                               f = f, &
                                               QLeft = f % storage(1) % Q(:,i,j), &
                                               QRight = f % storage(2) % Q(:,i,j), &
                                               U_xLeft = f % storage(1) % U_x(:,i,j), &
                                               U_yLeft = f % storage(1) % U_y(:,i,j), &
                                               U_zLeft = f % storage(1) % U_z(:,i,j), &
                                               U_xRight = f % storage(2) % U_x(:,i,j), &
                                               U_yRight = f % storage(2) % U_y(:,i,j), &
                                               U_zRight = f % storage(2) % U_z(:,i,j), &
                                               mu = mu, &
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

      call f % ProjectFluxToElements(NINC, fStar, (/1, HMESH_NONE/))

      END SUBROUTINE computeBoundaryFlux_iNS
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
         use PhysicsStorage, only: NINC
         implicit none
         type(HexMesh)                  :: mesh
         real(kind=RP),      intent(in) :: time
         procedure(BCState_FCN)         :: externalStateProcedure

         call ViscousDiscretization % ComputeGradient( NINC, NINC, mesh , time , externalStateProcedure, iNSGradientValuesForQ_0D, iNSGradientValuesForQ_3D)

      end subroutine DGSpatial_ComputeGradient

      subroutine DensityLimiter(N,Q)
         implicit none
         integer,       intent(in)    :: N(3)
         real(kind=RP), intent(inout) :: Q(1:NINC,0:N(1),0:N(2),0:N(3))
!
!        ---------------
!        Local variables
!        ---------------
!
         integer  :: i, j, k 
         real(kind=RP) :: rhoIn01, p, rhomin, rhomax

         rhomin = thermodynamics % rho_min / refValues % rho
         rhomax = thermodynamics % rho_max / refValues % rho

         do k = 0, N(3) ; do j = 0, N(2) ; do i = 0, N(1)
            rhoIn01 = (Q(INSRHO,i,j,k)-rhomin)/(rhomax-rhomin)

            if ( rhoIn01 .ge. 1.0_RP ) then
               Q(INSRHO,i,j,k) = rhomax

            elseif ( rhoIn01 .le. 0.0_RP ) then
               Q(INSRHO,i,j,k) = rhomin

            else
               !p = POW3(rhoIn01)*(6.0_RP*POW2(rhoIn01)-15.0_RP*rhoIn01+10.0_RP)
               p = rhoIn01

               Q(INSRHO,i,j,k) = (rhomax-rhomin)*p + rhomin
         
            end if

         end do         ; end do         ; end do

      end subroutine DensityLimiter
!
!////////////////////////////////////////////////////////////////////////////////////////
!
end module SpatialDiscretization
