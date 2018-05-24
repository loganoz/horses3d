!
!//////////////////////////////////////////////////////
!
!   @File:    SpatialDiscretization.f90
!   @Author:  Juan (juan.manzanero@upm.es)
!   @Created: Tue Apr 24 17:10:06 2018
!   @Last revision date: Thu May 24 12:03:15 2018
!   @Last revision author: Juan Manzanero (juan.manzanero@upm.es)
!   @Last revision commit: a9728294bcfa3ec9f4c553776074055792be41e2
!
!//////////////////////////////////////////////////////
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
      use FaceClass
      use ElementClass
      use PhysicsStorage
      use Physics
      use MPI_Face_Class
      use MPI_Process_Info
      use DGSEMClass
      use ParticlesClass
      use FluidData
      use VariableConversion
      use BoundaryConditionFunctions
      use GradientsStabilization
#ifdef _HAS_MPI_
      use mpi
#endif

      private
      public  ComputeLaplacian, DGSpatial_ComputeGradient
      public  Initialize_SpaceAndTimeMethods, ComputeTimeDerivative, ComputeTimeDerivativeIsolated
      public  ComputeTimeDerivative_onlyLinear, ComputetimeDerivative_onlyNonLinear
      public  Finalize_SpaceAndTimeMethods
      public  viscousDiscretizationKey, CHDiscretizationKey


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

      interface
         subroutine UserDefinedSourceTermNS(x, Q, time, S, thermodynamics_, dimensionless_, refValues_)
            use SMConstants
            USE HexMeshClass
            use PhysicsStorage
            use FluidData
            IMPLICIT NONE
            real(kind=RP),             intent(in) :: x(NDIM)
            real(kind=RP),             intent(in) :: Q(NCONS)
            REAL(KIND=RP),             intent(in) :: time
            real(kind=RP),             intent(in) :: S(NCONS)
            type(Thermodynamics_t),    intent(in) :: thermodynamics_
            type(Dimensionless_t),     intent(in) :: dimensionless_
            type(RefValues_t),         intent(in) :: refValues_
         end subroutine UserDefinedSourceTermNS
      end interface

      procedure(computeElementInterfaceFluxF), pointer :: computeElementInterfaceFlux => computeElementInterfaceFlux_NS
      procedure(computeMPIFaceFluxF),          pointer :: computeMPIFaceFlux          => computeMPIFaceFlux_NS
      procedure(computeBoundaryFluxF),         pointer :: computeBoundaryFlux         => computeBoundaryFlux_NS

      logical :: enable_speed = .true.

      character(len=LINE_LENGTH), parameter  :: viscousDiscretizationKey = "viscous discretization"
      character(len=LINE_LENGTH), parameter  :: CHDiscretizationKey      = "cahn-hilliard discretization"

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
         character(len=LINE_LENGTH)       :: CHDiscretizationName
         
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
            if ( flowIsNavierStokes ) then
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

               call ViscousDiscretization % Construct(controlVariables, ViscousFlux0D, ViscousFlux2D, ViscousFlux3D, GetNSCHViscosity)
               call ViscousDiscretization % Describe
      
            else
               if (.not. allocated(ViscousDiscretization)) allocate(EllipticDiscretization_t :: ViscousDiscretization)
               call ViscousDiscretization % Construct(controlVariables, ViscousFlux0D, ViscousFlux2D, ViscousFlux3D, GetNSCHViscosity)
               
            end if

!
!           Initialize models
!           -----------------
            call InitializeLESModel(LESModel, controlVariables)
         
!   
!           Initialize Cahn-Hilliard discretization
!           ---------------------------------------         
            if ( .not. controlVariables % ContainsKey(CHDiscretizationKey) ) then
               print*, "Input file is missing entry for keyword: viscous discretization"
               errorMessage(STD_OUT)
               stop
            end if

            CHDiscretizationName = controlVariables % stringValueForKey(CHDiscretizationKey, requestedLength = LINE_LENGTH)
            call toLower(CHDiscretizationName)
            
            select case ( trim(CHDiscretizationName) )
            case("br1")
               allocate(BassiRebay1_t     :: CHDiscretization)
   
            case("br2")
               allocate(BassiRebay2_t     :: CHDiscretization)
   
            case("ip")
               allocate(InteriorPenalty_t :: CHDiscretization)
   
            case default
               write(STD_OUT,'(A,A,A)') 'Requested viscous discretization "',trim(CHDiscretizationName),'" is not implemented.'
               write(STD_OUT,'(A)') "Implemented discretizations are:"
               write(STD_OUT,'(A)') "  * BR1"
               write(STD_OUT,'(A)') "  * BR2"
               write(STD_OUT,'(A)') "  * IP"
               errorMessage(STD_OUT)
               stop 
   
            end select
   
            call CHDiscretization % Construct(controlVariables, CHDivergenceFlux0D, CHDivergenceFlux2D, CHDivergenceFlux3D, GetCHViscosity)
            call CHDiscretization % Describe
         
         end if
!
!        Compute wall distances
!        ----------------------
        ! call mesh % ComputeWallDistances
         
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
         class(Element), pointer    :: e
         INTEGER :: k, eID, fID, i, j
!
!        *****************************
!        Obtain the NS time derivative
!        *****************************
!
!$omp single
         call mesh % SetStorageToEqn(1)
!$omp end single
!
!        -----------------------------------------
!        Prolongation of the solution to the faces
!        -----------------------------------------
!
!$omp parallel shared(mesh, time) private(k, eID, fID, i, j)
         call mesh % ProlongSolutionToFaces(NCONS)
!        ----------------
!        Update MPI Faces
!        ----------------
!
#ifdef _HAS_MPI_
!$omp single
errorMessage(STD_OUT)
stop
         !call mesh % UpdateMPIFacesSolution
!$omp end single
#endif
!
!        -----------------
!        Compute gradients
!        -----------------
!
         if ( computeGradients ) then
            call ViscousDiscretization % ComputeGradient( NCONS, NGRAD, mesh , time , BCFunctions(NS_BC) % externalState, NSGradientValuesForQ_0D, NSGradientValuesForQ_3D)
         end if

#ifdef _HAS_MPI_
!$omp single
         if ( flowIsNavierStokes ) then
errorMessage(STD_OUT)
stop
            !call mesh % UpdateMPIFacesGradients
         end if
!$omp end single
#endif
!
!        -----------------------
!        Compute time derivative
!        -----------------------
!
         call TimeDerivative_ComputeQDot(mesh              = mesh , &
                                         particles         = particles, &
                                         t                 = time, &
                                         externalState     = BCFunctions(NS_BC) % externalState, &
                                         externalGradients = BCFunctions(NS_BC) % externalGradients )
!
!        *****************************************
!        Compute the Cahn-Hilliard time derivative
!        *****************************************
!
!        ------------------------------
!        Change memory to concentration
!        ------------------------------
!
!$omp single
         call mesh % SetStorageToEqn(2)
!$omp end single
!
!        -----------------------------------------
!        Prolongation of the solution to the faces
!        -----------------------------------------
!
         call mesh % ProlongSolutionToFaces(NCOMP)
!
!        ----------------
!        Update MPI Faces
!        ----------------
!
#ifdef _HAS_MPI_
!$omp single
errorMessage(STD_OUT)
stop
         !call mesh % UpdateMPIFacesSolution
!$omp end single
#endif
!
!        -----------------
!        Compute gradients
!        -----------------
!
         call CHDiscretization % ComputeGradient( NCOMP, NCOMP, mesh , time , BCFunctions(C_BC) % externalState, CHGradientValuesForQ_0D, CHGradientValuesForQ_3D)

#ifdef _HAS_MPI_
!$omp single
errorMessage(STD_OUT)
stop
         !call mesh % UpdateMPIFacesGradients
!$omp end single
#endif
!
!        ------------------------------
!        Compute the chemical potential
!        ------------------------------
!
!        Linear part
!        -----------
         call ComputeLaplacian(mesh = mesh , &
                               t    = time, &
                  externalState     = BCFunctions(C_BC) % externalState, &
                  externalGradients = BCFunctions(C_BC) % externalGradients )

!$omp do schedule(runtime) private(e)
         do eID = 1, mesh % no_of_elements
            e => mesh % elements(eID)
            e % storage % mu = - POW2(multiphase % eps) * e % storage % QDot
            call AddQuarticDWPDerivative(e % storage % c, e % storage % mu)
         end do
!$omp end do

!$omp single
         call mesh % SetStorageToEqn(3)
!$omp end single
!
!        *************************
!        Compute cDot: Q stores mu
!        *************************
!
!        -----------------------------------------
!        Prolongation of the solution to the faces
!        -----------------------------------------
!
         call mesh % ProlongSolutionToFaces(NCOMP)
!
!        ----------------
!        Update MPI Faces
!        ----------------
!
#ifdef _HAS_MPI_
!$omp single
errorMessage(STD_OUT)
stop
         !call mesh % UpdateMPIFacesSolution
!$omp end single
#endif
!
!        -----------------
!        Compute gradients
!        -----------------
!
         call CHDiscretization % ComputeGradient( NCOMP, NCOMP, mesh , time , BCFunctions(MU_BC) % externalState, CHGradientValuesForQ_0D, CHGradientValuesForQ_3D)

#ifdef _HAS_MPI_
!$omp single
!         call mesh % UpdateMPIFacesGradients
errorMessage(STD_OUT)
stop
!$omp end single
#endif
!
!        ------------------------------
!        Compute the chemical potential
!        ------------------------------
!
         call ComputeLaplacian(mesh = mesh , &
                               t    = time, &
                  externalState     = BCFunctions(MU_BC) % externalState, &
                  externalGradients = BCFunctions(MU_BC) % externalGradients )
!
!        Scale QDot with the Peclet number
!        ---------------------------------
!$omp do schedule(runtime) private(e)
         do eID = 1, mesh % no_of_elements
            e => mesh % elements(eID)
            e % storage % QDot = (1.0_RP / multiphase % Pe) * e % storage % QDot
         end do
!$omp end do
!
!        *****************************
!        Return the concentration to Q
!        *****************************
!
!$omp single
         call mesh % SetStorageToEqn(2)
!$omp end single
!
!        ***********************************
!        Compute the concentration advection
!        ***********************************
!
         if ( enable_speed ) then
!
!        Perform the stabilization
!        -------------------------
         call StabilizeGradients(mesh, time, BCFunctions(C_BC) % externalState)
!
!        Add the velocity field
!        ----------------------
!$omp do schedule(runtime) private(e)
         do eID = 1, mesh % no_of_elements
            e => mesh % elements(eID)
         
            do k = 0, e % Nxyz(3) ; do j = 0, e % Nxyz(2)   ; do i = 0, e % Nxyz(1)
               e % storage % cDot(1,i,j,k) = e % storage % cDot(1,i,j,k) - (&
                                   e % storage % QNS(IRHOU,i,j,k) * e % storage % c_x(1,i,j,k) &
                                 + e % storage % QNS(IRHOV,i,j,k) * e % storage % c_y(1,i,j,k) &
                                 + e % storage % QNS(IRHOW,i,j,k) * e % storage % c_z(1,i,j,k) &
                                 ) / e % storage % QNS(IRHO,i,j,k) 
            end do                ; end do                  ; end do
         end do
!$omp end do
         end if
!
!        ****************************
!        Return NS as default storage
!        ****************************
!
!$omp single
         call mesh % SetStorageToEqn(1)
!$omp end single
!
!        ****************************
!        Compute the Capilar pressure
!        ****************************
!
!$omp do schedule(runtime) private(e)
         do eID = 1, mesh % no_of_elements
            e => mesh % elements(eID) 
            do k = 0, e % Nxyz(3) ; do j = 0, e % Nxyz(2)   ; do i = 0, e % Nxyz(1)
               e % storage % QDot(IRHOU,i,j,k) =   e % storage % QDot(IRHOU,i,j,k) &
                                                      + (1.0_RP / (dimensionless % Re * multiphase % Ca)) * e % storage % mu(1,i,j,k) * e % storage % c_x(1,i,j,k)
               e % storage % QDot(IRHOV,i,j,k) =   e % storage % QDot(IRHOV,i,j,k) &
                                                      + (1.0_RP / (dimensionless % Re * multiphase % Ca)) * e % storage % mu(1,i,j,k) * e % storage % c_y(1,i,j,k)
               e % storage % QDot(IRHOW,i,j,k) =   e % storage % QDot(IRHOW,i,j,k) &
                                                      + (1.0_RP / (dimensionless % Re * multiphase % Ca)) * e % storage % mu(1,i,j,k) * e % storage % c_z(1,i,j,k)

            end do                ; end do                  ; end do
         end do
!$omp end parallel
!
      END SUBROUTINE ComputeTimeDerivative

      subroutine ComputeTimeDerivative_OnlyLinear( mesh, particles, time, BCFunctions)
         IMPLICIT NONE 
!
!        ---------
!        Arguments
!        ---------
!
         TYPE(HexMesh), target      :: mesh
         type(Particles_t)               :: particles
         REAL(KIND=RP)              :: time
         type(BCFunctions_t), intent(in)  :: BCFunctions(no_of_BCsets)
!
!        ---------------
!        Local variables
!        ---------------
!
         INTEGER :: i, j, k, eID, fID
         class(Element), pointer  :: e
         class(Face),    pointer  :: f
!
!        **************************************
!        Compute chemical potential: Q stores c
!        **************************************
!
!$omp parallel shared(mesh, time) private(e, i, j, k, eID, fID)
!
!        ------------------------------
!        Change memory to concentration
!        ------------------------------
!
!$omp single
         call mesh % SetStorageToEqn(2)
!$omp end single
!
!        -----------------------------------------
!        Prolongation of the solution to the faces
!        -----------------------------------------
!
         call mesh % ProlongSolutionToFaces(NCOMP)
!
!        ----------------
!        Update MPI Faces
!        ----------------
!
#ifdef _HAS_MPI_
!$omp single
errorMessage(STD_OUT)
stop
         !call mesh % UpdateMPIFacesSolution
!$omp end single
#endif
!
!        -----------------
!        Compute gradients
!        -----------------
!
         call CHDiscretization % ComputeGradient( NCOMP, NCOMP, mesh , time , BCFunctions(C_BC) % externalState, CHGradientValuesForQ_0D, CHGradientValuesForQ_3D)

#ifdef _HAS_MPI_
!$omp single
errorMessage(STD_OUT)
stop
         !call mesh % UpdateMPIFacesGradients
!$omp end single
#endif
!
!        ------------------------------
!        Compute the chemical potential
!        ------------------------------
!
!        Linear part
!        -----------
         call ComputeLaplacian(mesh = mesh , &
                               t    = time, &
                  externalState     = BCFunctions(C_BC) % externalState, &
                  externalGradients = BCFunctions(C_BC) % externalGradients )

!$omp do schedule(runtime)
         do eID = 1, mesh % no_of_elements
            e => mesh % elements(eID)
            e % storage % mu = - POW2(multiphase % eps) * e % storage % QDot
         end do
!$omp end do

!$omp single
         call mesh % SetStorageToEqn(3)
!$omp end single
!
!        *************************
!        Compute cDot: Q stores mu
!        *************************
!
!        -----------------------------------------
!        Prolongation of the solution to the faces
!        -----------------------------------------
!
         call mesh % ProlongSolutionToFaces(NCOMP)
!
!        ----------------
!        Update MPI Faces
!        ----------------
!
#ifdef _HAS_MPI_
!$omp single
!         call mesh % UpdateMPIFacesSolution
errorMessage(STD_OUT)
stop
!$omp end single
#endif
!
!        -----------------
!        Compute gradients
!        -----------------
!
         call CHDiscretization % ComputeGradient(NCOMP, NCOMP, mesh , time , BCFunctions(MU_BC) % externalState, CHGradientValuesForQ_0D, CHGradientValuesForQ_3D)

#ifdef _HAS_MPI_
!$omp single
errorMessage(STD_OUT)
stop
         !call mesh % UpdateMPIFacesGradients
!$omp end single
#endif
!
!        ------------------------------
!        Compute the chemical potential
!        ------------------------------
!
 call ComputeLaplacian(mesh = mesh , &
                               t    = time, &
                  externalState     = BCFunctions(MU_BC) % externalState, &
                  externalGradients = BCFunctions(MU_BC) % externalGradients )
!
!        Scale QDot with the Peclet number
!        ---------------------------------
!$omp do schedule(runtime)
         do eID = 1, mesh % no_of_elements
            e => mesh % elements(eID)
            e % storage % QDot = (1.0_RP / multiphase % Pe) * e % storage % QDot
         end do
!$omp end do
!
!        *****************************
!        Return the concentration to Q
!        *****************************
!
!$omp single
         call mesh % SetStorageToEqn(2)
!$omp end single
!$omp end parallel

      end subroutine ComputeTimeDerivative_OnlyLinear

      subroutine ComputeTimeDerivative_OnlyNonLinear( mesh, particles, time, BCFunctions)
         IMPLICIT NONE 
!
!        ---------
!        Arguments
!        ---------
!
         TYPE(HexMesh), target      :: mesh
         type(Particles_t)               :: particles
         REAL(KIND=RP)              :: time
         type(BCFunctions_t), intent(in)  :: BCFunctions(no_of_BCsets)
!
!        ---------------
!        Local variables
!        ---------------
!
         INTEGER :: i, j, k, eID, fID
         class(Element), pointer  :: e
         class(Face),    pointer  :: f

!
!        **************************************
!        Compute chemical potential: Q stores c
!        **************************************
!
!$omp parallel shared(mesh, time) private(e, i, j, k, eID, fID)
!
!        ------------------------------
!        Change memory to concentration
!        ------------------------------
!
!$omp single
         call mesh % SetStorageToEqn(2)
!$omp end single
!
!        -----------------------------------------
!        Prolongation of the solution to the faces
!        -----------------------------------------
!
         call mesh % ProlongSolutionToFaces(NCOMP)
!
!        ----------------
!        Update MPI Faces
!        ----------------
!
#ifdef _HAS_MPI_
!$omp single
errorMessage(STD_OUT)
stop
         !call mesh % UpdateMPIFacesSolution
!$omp end single
#endif
!
!        -----------------
!        Compute gradients
!        -----------------
!
         call CHDiscretization % ComputeGradient( NCOMP, NCOMP, mesh , time , BCFunctions(C_BC) % externalState, CHGradientValuesForQ_0D, CHGradientValuesForQ_3D)

#ifdef _HAS_MPI_
!$omp single
errorMessage(STD_OUT)
stop
         !call mesh % UpdateMPIFacesGradients
!$omp end single
#endif
!
!        ------------------------------
!        Compute the chemical potential
!        ------------------------------
!
!        Linear part: only Neumann boundary conditions contribution
!        -----------
         call ComputeLaplacianNeumannBCs(mesh = mesh , &
                               t    = time, &
                  externalState     = BCFunctions(C_BC) % externalState, &
                  externalGradients = BCFunctions(C_BC) % externalGradients )

!$omp do schedule(runtime) private(e)
         do eID = 1, mesh % no_of_elements
            e => mesh % elements(eID)
            e % storage % mu = - POW2(multiphase % eps) * e % storage % QDot
            call AddQuarticDWPDerivative(e % storage % c, e % storage % mu)
         end do
!$omp end do

!$omp single
         call mesh % SetStorageToEqn(3)
!$omp end single
!
!        *************************
!        Compute cDot: Q stores mu
!        *************************
!
!        -----------------------------------------
!        Prolongation of the solution to the faces
!        -----------------------------------------
!
         call mesh % ProlongSolutionToFaces(NCOMP)
!
!        ----------------
!        Update MPI Faces
!        ----------------
!
#ifdef _HAS_MPI_
!$omp single
errorMessage(STD_OUT)
stop
         !call mesh % UpdateMPIFacesSolution
!$omp end single
#endif
!
!        -----------------
!        Compute gradients
!        -----------------
!
         call CHDiscretization % ComputeGradient( NCOMP, NCOMP, mesh , time , BCFunctions(MU_BC) % externalState, CHGradientValuesForQ_0D, CHGradientValuesForQ_3D)

#ifdef _HAS_MPI_
!$omp single
!         call mesh % UpdateMPIFacesGradients
errorMessage(STD_OUT)
stop
!$omp end single
#endif
!
!        ------------------------------
!        Compute the chemical potential
!        ------------------------------
!
         call ComputeLaplacian(mesh = mesh , &
                               t    = time, &
                  externalState     = BCFunctions(MU_BC) % externalState, &
                  externalGradients = BCFunctions(MU_BC) % externalGradients )
!
!        Scale QDot with the Peclet number
!        ---------------------------------
!$omp do schedule(runtime) private(e)
         do eID = 1, mesh % no_of_elements
            e => mesh % elements(eID)
            e % storage % QDot = (1.0_RP / multiphase % Pe) * e % storage % QDot
         end do
!$omp end do
!
!        *****************************
!        Return the concentration to Q
!        *****************************
!
!$omp single
         call mesh % SetStorageToEqn(2)
!$omp end single
!
!        ***********************************
!        Compute the concentration advection
!        ***********************************
!
         if ( enable_speed ) then
!
!        Perform the stabilization
!        -------------------------
         call StabilizeGradients(mesh, time, BCFunctions(C_BC) % externalState)
!
!        Add the velocity field
!        ----------------------
!$omp do schedule(runtime)
         do eID = 1, mesh % no_of_elements
            e => mesh % elements(eID)
         
            do k = 0, e % Nxyz(3) ; do j = 0, e % Nxyz(2)   ; do i = 0, e % Nxyz(1)
               e % storage % cDot(1,i,j,k) = e % storage % cDot(1,i,j,k) - (&
                                   e % storage % QNS(IRHOU,i,j,k) * e % storage % c_x(1,i,j,k) &
                                 + e % storage % QNS(IRHOV,i,j,k) * e % storage % c_y(1,i,j,k) &
                                 + e % storage % QNS(IRHOW,i,j,k) * e % storage % c_z(1,i,j,k) &
                                 ) / e % storage % QNS(IRHO,i,j,k) 
            end do                ; end do                  ; end do
         end do
!$omp end do
         end if
!$omp end parallel

      end subroutine ComputeTimeDerivative_OnlyNonLinear
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
!$omp do schedule(runtime) private(eID,fID,i,j,k)
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
errorMessage(STD_OUT)
stop
            if ( flowIsNavierStokes ) then 
!               call mesh % GatherMPIFacesGradients
            else  
!               call mesh % GatherMPIFacesSolution
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
         if (.not. mesh % child) then
!$omp do schedule(runtime) private(i,j,k)
            do eID = 1, mesh % no_of_elements
               associate ( e => mesh % elements(eID) )
               do k = 0, e % Nxyz(3)   ; do j = 0, e % Nxyz(2) ; do i = 0, e % Nxyz(1)
                  call UserDefinedSourceTermNS(e % geom % x(:,i,j,k), e % storage % Q(:,i,j,k), t, e % storage % S_NS(:,i,j,k), thermodynamics, dimensionless, refValues)
               end do                  ; end do                ; end do
               end associate
            end do
!$omp end do
         end if

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
               e % storage % QDot(:,i,j,k) = e % storage % QDot(:,i,j,k) + e % storage % S_NS(:,i,j,k)
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
            call ViscousDiscretization  % ComputeInnerFluxes ( NCONS, NGRAD, e , viscousContravariantFlux) 

         else
!
!           With LES model
!           --------------
            call ViscousDiscretization  % ComputeInnerFluxesWithSGS ( e , viscousContravariantFlux  ) 

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
            contravariantFlux = inviscidContravariantFlux - viscousContravariantFlux
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
         real(kind=RP) :: mu, muL, muR

         if ( .not. LESModel % active ) then
         DO j = 0, f % Nf(2)
            DO i = 0, f % Nf(1)

!      
!              --------------
!              Viscous fluxes
!              --------------
!      
               call ViscousDiscretization % GetViscosity(f % storage(1) % c(1,i,j), muL)
               call ViscousDiscretization % GetViscosity(f % storage(2) % c(1,i,j), muR)
               mu = 0.5_RP * (muL + muR)

               CALL ViscousDiscretization % RiemannSolver(nEqn = NCONS, nGradEqn = NGRAD, &
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
               CALL ViscousDiscretization % RiemannSolverWithSGS(f = f, &
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
         real(kind=RP) :: mu, muL, muR
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
               call ViscousDiscretization % GetViscosity(f % storage(1) % c(1,i,j), muL)
               call ViscousDiscretization % GetViscosity(f % storage(2) % c(1,i,j), muR)
               mu = 0.5_RP * (muL + muR)

               CALL ViscousDiscretization % RiemannSolver(nEqn = NCONS, nGradEqn = NGRAD, &
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
      real(kind=RP)                   :: mu, muL, muR
      
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
            call ViscousDiscretization % GetViscosity(f % storage(1) % c(1,i,j), muL)
            call ViscousDiscretization % GetViscosity(f % storage(2) % c(1,i,j), muR)
            mu = 0.5_RP * (muL + muR)

            CALL ViscousDiscretization % RiemannSolver(nEqn = NCONS, nGradEqn = NGRAD, &
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

      call f % ProjectFluxToElements(NCONS, fStar, (/1, HMESH_NONE/))

      END SUBROUTINE computeBoundaryFlux_NS
!
!///////////////////////////////////////////////////////////////////////////////////////////
!
!           Procedures to compute the state variables Laplacian
!           ---------------------------------------------------
!
!///////////////////////////////////////////////////////////////////////////////////////////
!
      subroutine ComputeLaplacian( mesh , t, externalState, externalGradients )
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
!
!        ****************
!        Volume integrals
!        ****************
!
!$omp do schedule(runtime) 
         do eID = 1 , size(mesh % elements)
            call Laplacian_VolumetricContribution( mesh % elements(eID) , t)
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
               CALL computeElementInterfaceFlux_Laplacian( f ) 
 
            case (HMESH_BOUNDARY) 
               CALL computeBoundaryFlux_Laplacian(f, t, externalState, externalGradients) 
            
            case (HMESH_MPI)

            case default
               print*, "Unrecognized face type"
               errorMessage(STD_OUT)
               stop
                
            end select 
            end associate 
         end do 
!$omp end do 
!
!        ***************************************************************
!        Surface integrals and scaling of elements with non-shared faces
!        ***************************************************************
! 
!$omp do schedule(runtime) private(i, j, k)
         do eID = 1, size(mesh % elements) 
            associate(e => mesh % elements(eID)) 
            if ( e % hasSharedFaces ) cycle
            call Laplacian_FacesContribution(e, t, mesh) 
 
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
errorMessage(STD_OUT)
stop
!            call mesh % GatherMPIFacesGradients
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
                  CALL computeMPIFaceFlux_Laplacian( f ) 
               end select 
               end associate 
            end do 
!$omp end do 
!
!           ***********************************************************
!           Surface integrals and scaling of elements with shared faces
!           ***********************************************************
! 
!$omp do schedule(runtime) 
            do eID = 1, size(mesh % elements) 
               associate(e => mesh % elements(eID)) 
               if ( .not. e % hasSharedFaces ) cycle
               call Laplacian_FacesContribution(e, t, mesh) 
 
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

      end subroutine ComputeLaplacian

      subroutine ComputeLaplacianNeumannBCs( mesh , t, externalState, externalGradients )
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
!
!        **************************
!        Reset QDot and face fluxes
!        **************************
!
         do eID = 1, mesh % no_of_elements
            mesh % elements(eID) % storage % QDot = 0.0_RP
         end do
   
         do fID = 1, size(mesh % faces)
            mesh % faces(fID) % storage(1) % genericInterfaceFluxMemory = 0.0_RP
            mesh % faces(fID) % storage(2) % genericInterfaceFluxMemory = 0.0_RP
         end do
!
!        ******************************************
!        Compute Riemann solver of non-shared faces
!        ******************************************
!
!$omp do schedule(runtime) 
         do fID = 1, size(mesh % faces) 
            associate( f => mesh % faces(fID)) 
            select case (f % faceType) 
            case (HMESH_BOUNDARY) 
               CALL computeBoundaryFlux_Laplacian(f, t, externalState, externalGradients) 
            end select 
            end associate 
         end do 
!$omp end do 
!
!        ***************************************************************
!        Surface integrals and scaling of elements with non-shared faces
!        ***************************************************************
! 
!$omp do schedule(runtime) private(i, j, k)
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
            !call mesh % GatherMPIFacesGradients
errorMessage(STD_OUT)
stop
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
!           Add a MPI Barrier
!           -----------------
!$omp single
            call mpi_barrier(MPI_COMM_WORLD, ierr)
!$omp end single
         end if
#endif

      end subroutine ComputeLaplacianNeumannBCs

!
!///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
!
      subroutine Laplacian_VolumetricContribution( e , t )
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
         real(kind=RP) :: contravariantFlux  ( 1:NCOMP, 0:e%Nxyz(1) , 0:e%Nxyz(2) , 0:e%Nxyz(3), 1:NDIM ) 
         integer       :: eID
!
!        *************************************
!        Compute interior contravariant fluxes
!        *************************************
!
!        Compute contravariant flux
!        --------------------------
         call CHDiscretization  % ComputeInnerFluxes (NCOMP, NCOMP, e , contravariantFlux  ) 
!
!        ************************
!        Perform volume integrals
!        ************************
!
         e % storage % QDot = - ScalarWeakIntegrals % StdVolumeGreen ( e , NCOMP, contravariantFlux ) 

      end subroutine Laplacian_VolumetricContribution
!
!///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
!
      subroutine Laplacian_FacesContribution( e , t , mesh)
         use HexMeshClass
         use PhysicsStorage
         implicit none
         type(Element)           :: e
         real(kind=RP)           :: t
         type(HexMesh)           :: mesh

         e % storage % QDot = e % storage % QDot + ScalarWeakIntegrals % StdFace(e, NCOMP, &
                      mesh % faces(e % faceIDs(EFRONT))  % storage(e % faceSide(EFRONT))  % fStar, &
                      mesh % faces(e % faceIDs(EBACK))   % storage(e % faceSide(EBACK))   % fStar, &
                      mesh % faces(e % faceIDs(EBOTTOM)) % storage(e % faceSide(EBOTTOM)) % fStar, &
                      mesh % faces(e % faceIDs(ERIGHT))  % storage(e % faceSide(ERIGHT))  % fStar, &
                      mesh % faces(e % faceIDs(ETOP))    % storage(e % faceSide(ETOP))    % fStar, &
                      mesh % faces(e % faceIDs(ELEFT))   % storage(e % faceSide(ELEFT))   % fStar )

      end subroutine Laplacian_FacesContribution
!
!///////////////////////////////////////////////////////////////////////////////////////////// 
! 
!        Riemann solver drivers 
!        ---------------------- 
! 
!///////////////////////////////////////////////////////////////////////////////////////////// 
! 
      subroutine computeElementInterfaceFlux_Laplacian(f)
         use FaceClass
         use Physics
         use PhysicsStorage
         IMPLICIT NONE
         TYPE(Face)   , INTENT(inout) :: f   
         integer       :: i, j
         real(kind=RP) :: flux(1:NCOMP,0:f % Nf(1),0:f % Nf(2))

         DO j = 0, f % Nf(2)
            DO i = 0, f % Nf(1)

!      
!              --------------
!              Viscous fluxes
!              --------------
!      
               CALL CHDiscretization % RiemannSolver(nEqn = NCOMP, nGradEqn = NCOMP, &
                                                  f = f, &
                                                  QLeft = f % storage(1) % Q(:,i,j), &
                                                  QRight = f % storage(2) % Q(:,i,j), &
                                                  U_xLeft = f % storage(1) % U_x(:,i,j), &
                                                  U_yLeft = f % storage(1) % U_y(:,i,j), &
                                                  U_zLeft = f % storage(1) % U_z(:,i,j), &
                                                  U_xRight = f % storage(2) % U_x(:,i,j), &
                                                  U_yRight = f % storage(2) % U_y(:,i,j), &
                                                  U_zRight = f % storage(2) % U_z(:,i,j), &
                                                  mu = 1.0_RP, &
                                                  nHat = f % geom % normal(:,i,j) , &
                                                  dWall = f % geom % dWall(i,j), &
                                                  flux  = flux(:,i,j) )

               flux(:,i,j) = flux(:,i,j) * f % geom % jacobian(i,j)

            end do
         end do
!
!        ---------------------------
!        Return the flux to elements
!        ---------------------------
!
         call f % ProjectFluxToElements(NCOMP, flux, (/1,2/))

      end subroutine computeElementInterfaceFlux_Laplacian

      subroutine computeMPIFaceFlux_Laplacian(f)
         use FaceClass
         use Physics
         use PhysicsStorage
         IMPLICIT NONE
         TYPE(Face)   , INTENT(inout) :: f   
         integer       :: i, j
         integer       :: thisSide
         real(kind=RP) :: flux(1:NCOMP,0:f % Nf(1),0:f % Nf(2))

         DO j = 0, f % Nf(2)
            DO i = 0, f % Nf(1)
!      
!              --------------
!              Viscous fluxes
!              --------------
!      
               CALL CHDiscretization % RiemannSolver(nEqn = NCOMP, nGradEqn = NCOMP, &
                                                  f = f, &
                                                  QLeft = f % storage(1) % Q(:,i,j), &
                                                  QRight = f % storage(2) % Q(:,i,j), &
                                                  U_xLeft = f % storage(1) % U_x(:,i,j), &
                                                  U_yLeft = f % storage(1) % U_y(:,i,j), &
                                                  U_zLeft = f % storage(1) % U_z(:,i,j), &
                                                  U_xRight = f % storage(2) % U_x(:,i,j), &
                                                  U_yRight = f % storage(2) % U_y(:,i,j), &
                                                  U_zRight = f % storage(2) % U_z(:,i,j), &
                                                  mu = 1.0_RP, &
                                                  nHat = f % geom % normal(:,i,j) , &
                                                  dWall = f % geom % dWall(i,j), &
                                                  flux  = flux(:,i,j) )

               flux(:,i,j) = flux(:,i,j) * f % geom % jacobian(i,j)

            END DO   
         END DO  
!
!        ---------------------------
!        Return the flux to elements: The sign in eR % storage % FstarB has already been accouted.
!        ---------------------------
!
         thisSide = maxloc(f % elementIDs, dim = 1)
         call f % ProjectFluxToElements(NCOMP, flux, (/thisSide, HMESH_NONE/))

      end subroutine ComputeMPIFaceFlux_Laplacian

      subroutine computeBoundaryFlux_Laplacian(f, time, externalStateProcedure , externalGradientsProcedure)
      USE ElementClass
      use FaceClass
      USE EllipticDiscretizations
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
      procedure(BCState_FCN)     :: externalState
      procedure(BCGradients_FCN) :: externalGradients
!
!     ---------------
!     Local variables
!     ---------------
!
      INTEGER                         :: i, j
      INTEGER, DIMENSION(2)           :: N
      REAL(KIND=RP)                   :: UGradExt(NDIM , NCOMP) 
      real(kind=RP)                   :: flux(NCOMP, 0:f % Nf(1), 0:f % Nf(2))
      
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

      do j = 0, f % Nf(2)  ; do i = 0, f % Nf(1)
         UGradExt(IX,:) = f % storage(1) % U_x(:,i,j)
         UGradExt(IY,:) = f % storage(1) % U_y(:,i,j)
         UGradExt(IZ,:) = f % storage(1) % U_z(:,i,j)

         CALL externalGradientsProcedure(  f % geom % x(:,i,j), &
                                           time, &
                                           f % geom % normal(:,i,j), &
                                           UGradExt,&
                                           boundaryType )

         f % storage(1) % U_x(:,i,j) = UGradExt(IX,:)
         f % storage(1) % U_y(:,i,j) = UGradExt(IY,:)
         f % storage(1) % U_z(:,i,j) = UGradExt(IZ,:)

         f % storage(2) % U_x(:,i,j) = UGradExt(IX,:)
         f % storage(2) % U_y(:,i,j) = UGradExt(IY,:)
         f % storage(2) % U_z(:,i,j) = UGradExt(IZ,:)

!   
!           --------------
!           Viscous fluxes
!           --------------
!   
         CALL CHDiscretization % RiemannSolver(nEqn = NCOMP, nGradEqn = NCOMP, &
                                            f = f, &
                                            QLeft = f % storage(1) % Q(:,i,j), &
                                            QRight = f % storage(2) % Q(:,i,j), &
                                            U_xLeft = f % storage(1) % U_x(:,i,j), &
                                            U_yLeft = f % storage(1) % U_y(:,i,j), &
                                            U_zLeft = f % storage(1) % U_z(:,i,j), &
                                            U_xRight = f % storage(2) % U_x(:,i,j), &
                                            U_yRight = f % storage(2) % U_y(:,i,j), &
                                            U_zRight = f % storage(2) % U_z(:,i,j), &
                                            mu = 1.0_RP, &
                                            nHat = f % geom % normal(:,i,j) , &
                                            dWall = f % geom % dWall(i,j), &
                                            flux  = flux(:,i,j) )

         flux(:,i,j) = flux(:,i,j) * f % geom % jacobian(i,j)

      end do               ; end do

      call f % ProjectFluxToElements(NCOMP, flux, (/1, HMESH_NONE/))

      end subroutine computeBoundaryFlux_Laplacian
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
         TYPE(Face)   , INTENT(inout) :: f
         integer      , intent(in)    :: thisSide
      END SUBROUTINE computeIsolatedFaceFlux_NS
!
!////////////////////////////////////////////////////////////////////////////////////////
!
end module SpatialDiscretization
