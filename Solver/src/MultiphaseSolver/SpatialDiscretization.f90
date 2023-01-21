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
      use FluidData
      use VariableConversion, only: mGradientVariables, GetmOneFluidViscosity,&
                                    GetmTwoFluidsViscosity, chGradientVariables,&
                                    GetCHViscosity
      use BoundaryConditions, only: BCs, SetBoundaryConditionsEqn, NS_BC, C_BC, MU_BC
      use ProblemFileFunctions, only: UserDefinedSourceTermNS_f
      use ParticlesClass
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

         SUBROUTINE computeBoundaryFluxF(f, time)
            use SMConstants
            use FaceClass,  only: Face
            IMPLICIT NONE
            type(Face),    intent(inout) :: f
            REAL(KIND=RP)                :: time
         end subroutine computeBoundaryFluxF
      end interface
      
      character(len=LINE_LENGTH), parameter  :: viscousDiscretizationKey = "viscous discretization"
      character(len=LINE_LENGTH), parameter     :: CHDiscretizationKey = "cahn-hilliard discretization"

      real(kind=RP), protected :: IMEX_S0 = 0.0_RP 
      real(kind=RP), protected :: IMEX_K0 = 1.0_RP
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
         character(len=LINE_LENGTH) :: inviscidDiscretizationName
         character(len=LINE_LENGTH) :: viscousDiscretizationName
         character(len=LINE_LENGTH) :: CHDiscretizationName
         
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

            case default
               write(STD_OUT,'(A,A,A)') 'Requested inviscid discretization "',trim(inviscidDiscretizationName),'" is not implemented.'
               write(STD_OUT,'(A)') "Implemented discretizations are:"
               write(STD_OUT,'(A)') "  * Standard"
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

               call ViscousDiscretization % Construct(controlVariables, ELLIPTIC_MU)
               call ViscousDiscretization % Describe
!
!           Compute wall distances
!           ----------------------
            call mesh % ComputeWallDistances

!
!           Initialize Cahn--Hilliard discretization
!           ----------------------------------------         
            if ( .not. controlVariables % ContainsKey(CHDiscretizationKey) ) then
               print*, "Input file is missing entry for keyword: Cahn-Hilliard discretization"
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
   
            call CHDiscretization % Construct(controlVariables, ELLIPTIC_CH)
            call CHDiscretization % Describe

         
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
      SUBROUTINE ComputeTimeDerivative( mesh, particles, time, mode)
         IMPLICIT NONE 
!
!        ---------
!        Arguments
!        ---------
!
         TYPE(HexMesh), target           :: mesh
         type(Particles_t)               :: particles
         REAL(KIND=RP)                   :: time
         integer, intent(in)             :: mode
!
!        ---------------
!        Local variables
!        ---------------
!
         INTEGER                 :: k, eID, fID, i, j
         real(kind=RP)           :: sqrtRho
         class(Element), pointer :: e


!$omp parallel shared(mesh, time)
!
!///////////////////////////////////////////////////
!        1st step: Get chemical potential
!///////////////////////////////////////////////////
!
!        ------------------------------------------
!        Update concentration with the state vector
!        ------------------------------------------
!
         select case (mode)
         case (CTD_IGNORE_MODE,CTD_IMEX_EXPLICIT)
!$omp do schedule(runtime)
            do eID = 1, size(mesh % elements)
               mesh % elements(eID) % storage % c(1,:,:,:) = mesh % elements(eID) % storage % QNS(IMC,:,:,:)
            end do
!$omp end do         
         end select

!
!        -------------------------------
!        Set memory to Cahn-Hilliard (C)
!        -------------------------------
!
!$omp single
         call mesh % SetStorageToEqn(C_BC)
         select case (mode)
         case (CTD_IGNORE_MODE,CTD_IMEX_EXPLICIT)
            call SetBoundaryConditionsEqn(C_BC)

         case (CTD_IMEX_IMPLICIT,CTD_LAPLACIAN)
            call SetBoundaryConditionsEqn(MU_BC)

         end select
!$omp end single
!
!        --------------------------------------------
!        Prolong Cahn-Hilliard concentration to faces
!        --------------------------------------------
!
         call mesh % ProlongSolutionToFaces(NCOMP)
!
!        ----------------
!        Update MPI Faces
!        ----------------
!
#ifdef _HAS_MPI_
!$omp single
         call mesh % UpdateMPIFacesSolution(NCOMP)
!$omp end single
#endif
!
!        ------------------------------------------------------------
!        Get concentration (lifted) gradients (also prolong to faces)
!        ------------------------------------------------------------
!
         call CHDiscretization % ComputeGradient(NCOMP, NCOMP, mesh, time, chGradientVariables)
!
!        --------------------
!        Update MPI Gradients
!        --------------------
!
#ifdef _HAS_MPI_
!$omp single
         call mesh % UpdateMPIFacesGradients(NCOMP)
!$omp end single
#endif
!
!        ----------------------
!        Get chemical potential
!        ----------------------
!
!        Get the concentration Laplacian (into QDot => cDot)

         call ComputeLaplacian(mesh, time)

         select case (mode)
         case (CTD_IGNORE_MODE, CTD_IMEX_EXPLICIT)
!$omp do schedule(runtime)
            do eID = 1, size(mesh % elements)
!
!            + Linear part
               !mesh % elements(eID) % storage % mu = - POW2(multiphase % eps)* mesh % elements(eID) % storage % QDot
               mesh % elements(eID) % storage % mu = - 1.5_RP * multiphase % eps * multiphase % sigma * mesh % elements(eID) % storage % QDot
!
!            + NonLinear part
               !call AddQuarticDWPDerivative(mesh % elements(eID) % storage % c, mesh % elements(eID) % storage % mu)
               call Multiphase_AddChemFEDerivative(mesh % elements(eID) % storage % c, mesh % elements(eID) % storage % mu)

            end do
!$omp end do         
         case (CTD_IMEX_IMPLICIT)
!$omp do schedule(runtime)
            do eID = 1, size(mesh % elements)
!
!            + Linear part
               !mesh % elements(eID) % storage % mu = - IMEX_K0 * POW2(multiphase % eps) * mesh % elements(eID) % storage % QDot &
               mesh % elements(eID) % storage % mu = - 1.5_RP * IMEX_K0 * multiphase % eps * multiphase % sigma * mesh % elements(eID) % storage % QDot &
                                                     + IMEX_S0 * mesh % elements(eID) % storage % c
!            + Multiply by mobility
               mesh % elements(eID) % storage % mu = multiphase % M0 * mesh % elements(eID) % storage % mu
            end do
!$omp end do         
         end select
!
!        -----------------------------------
!        Prolong chemical potential to faces
!        -----------------------------------
!
         select case(mode)
         case(CTD_LAPLACIAN)
         case default
!$omp single
         call mesh % SetStorageToEqn(MU_BC)
!$omp end single
         call mesh % ProlongSolutionToFaces(NCOMP)
!
!        ----------------
!        Update MPI Faces
!        ----------------
!
#ifdef _HAS_MPI_
!$omp single
         call mesh % UpdateMPIFacesSolution(NCOMP)
!$omp end single
#endif
         end select
!
!/////////////////////////////////////////////////////////////////////////////////
!        2nd step: If IMEX_IMPLCIIT, get the chemical potential laplacian and exit
!/////////////////////////////////////////////////////////////////////////////////
!
         select case (mode)
         case (CTD_IMEX_IMPLICIT)
!
!           ------------------------------------------------------------
!           Get concentration (lifted) gradients (also prolong to faces)
!           ------------------------------------------------------------
!
            call CHDiscretization % ComputeGradient(NCOMP, NCOMP, mesh, time, chGradientVariables)
!
!           --------------------
!           Update MPI Gradients
!           --------------------
!
#ifdef _HAS_MPI_
!$omp single
            call mesh % UpdateMPIFacesGradients(NCOMP)
!$omp end single
#endif
!
!           ----------------------
!           Get chemical potential
!           ----------------------
!
!           Get the concentration Laplacian (into QDot => cDot)

            call ComputeLaplacian(mesh, time)
!
!           ------------------------------------------
!           *** WARNING! The storage leaves set to CH!
!           ------------------------------------------
!
!$omp single
            call mesh % SetStorageToEqn(C_BC)
            call SetBoundaryConditionsEqn(C_BC)
!$omp end single
         end select
!
!///////////////////////////////////////////////
!        3rd step: Navier-Stokes time derivative
!///////////////////////////////////////////////
!
         select case (mode)
         case (CTD_IGNORE_MODE, CTD_IMEX_EXPLICIT)
!$omp single         
         call mesh % SetStorageToEqn(NS_BC)
         call SetBoundaryConditionsEqn(NS_BC)
!$omp end single
!
!        -------------------------
!        Prolong solution to faces        
!        -------------------------
!
         call mesh % ProlongSolutionToFaces(NCONS)
!
!        ----------------
!        Update MPI Faces
!        ----------------
!
#ifdef _HAS_MPI_
!$omp single
         call mesh % UpdateMPIFacesSolution(NCONS)
!$omp end single
#endif
!
!        -------------------------------------
!        Get the density in faces and elements
!        -------------------------------------
!
!$omp do schedule(runtime)
         do eID = 1, size(mesh % elements)
            mesh % elements(eID) % storage % rho = dimensionless % rho(2) + (dimensionless % rho(1)-dimensionless % rho(2))*mesh % elements(eID) % storage % Q(IMC,:,:,:)

            mesh % elements(eID) % storage % rho = min(max(mesh % elements(eID) % storage % rho, dimensionless % rho_min),dimensionless % rho_max)
         end do
!$omp end do nowait

!$omp do schedule(runtime)
         do fID = 1, size(mesh % faces)
            mesh % faces(fID) % storage(1) % rho = dimensionless % rho(2) + (dimensionless % rho(1)-dimensionless % rho(2))*mesh % faces(fID) % storage(1) % Q(IMC,:,:)
            mesh % faces(fID) % storage(2) % rho = dimensionless % rho(2) + (dimensionless % rho(1)-dimensionless % rho(2))*mesh % faces(fID) % storage(2) % Q(IMC,:,:)

            mesh % faces(fID) % storage(1) % rho = min(max(mesh % faces(fID) % storage(1) % rho, dimensionless % rho_min),dimensionless % rho_max)
            mesh % faces(fID) % storage(2) % rho = min(max(mesh % faces(fID) % storage(2) % rho, dimensionless % rho_min),dimensionless % rho_max)
         end do
!$omp end do
!
!        ----------------------------------------
!        Compute local entropy variables gradient
!        ----------------------------------------
!
         call ViscousDiscretization % ComputeLocalGradients( NCONS, NCONS, mesh , time , mGradientVariables)
!
!        --------------------
!        Update MPI Gradients
!        --------------------
!
#ifdef _HAS_MPI_
!$omp single
         call mesh % UpdateMPIFacesGradients(NCONS)
!$omp end single
#endif
!
!        -------------------------------------
!        Add the Non-Conservative term to QDot
!        -------------------------------------
!
!$omp do schedule(runtime) private(i,j,k,e,sqrtRho)
         do eID = 1, size(mesh % elements)
            e => mesh % elements(eID)
            do k = 0, e % Nxyz(3) ; do j = 0, e % Nxyz(2) ; do i = 0, e % Nxyz(1)
               sqrtRho = sqrt(e % storage % rho(i,j,k))
               e % storage % QDot(IMC,i,j,k)      = 0.0_RP
               e % storage % QDot(IMSQRHOU,i,j,k) = -0.5_RP*sqrtRho*(  e % storage % Q(IMSQRHOU,i,j,k)*e % storage % U_x(IGU,i,j,k) & 
                                                                     + e % storage % Q(IMSQRHOV,i,j,k)*e % storage % U_y(IGU,i,j,k) &   
                                                                     + e % storage % Q(IMSQRHOW,i,j,k)*e % storage % U_z(IGU,i,j,k) ) &
                                                    - e % storage % Q(IMC,i,j,k)*e % storage % U_x(IGMU,i,j,k)

               e % storage % QDot(IMSQRHOV,i,j,k) = -0.5_RP*sqrtRho*(  e % storage % Q(IMSQRHOU,i,j,k)*e % storage % U_x(IGV,i,j,k) & 
                                                                     + e % storage % Q(IMSQRHOV,i,j,k)*e % storage % U_y(IGV,i,j,k) &   
                                                                     + e % storage % Q(IMSQRHOW,i,j,k)*e % storage % U_z(IGV,i,j,k) ) &
                                                    - e % storage % Q(IMC,i,j,k)*e % storage % U_y(IGMU,i,j,k)

               e % storage % QDot(IMSQRHOW,i,j,k) = -0.5_RP*sqrtRho*(  e % storage % Q(IMSQRHOU,i,j,k)*e % storage % U_x(IGW,i,j,k) & 
                                                                     + e % storage % Q(IMSQRHOV,i,j,k)*e % storage % U_y(IGW,i,j,k) &   
                                                                     + e % storage % Q(IMSQRHOW,i,j,k)*e % storage % U_z(IGW,i,j,k) ) &
                                                    - e % storage % Q(IMC,i,j,k)*e % storage % U_z(IGMU,i,j,k)
   
               e % storage % QDot(IMP,i,j,k) = -dimensionless % invMa2*(  e % storage % U_x(IGU,i,j,k) + e % storage % U_y(IGV,i,j,k) &
                                                                          + e % storage % U_z(IGW,i,j,k))  

               e % storage % QDot(:,i,j,k) = e % storage % QDot(:,i,j,k) * e % geom % jacobian(i,j,k)
            end do                ; end do                ; end do

         end do
!$omp end do

         call ViscousDiscretization % LiftGradients( NCONS, NCONS, mesh , time , mGradientVariables)
!
!        -----------------------
!        Compute time derivative
!        -----------------------
!
         select case (mode)
         case(CTD_IMEX_EXPLICIT)
            call multiphase % SetStarMobility(0.0_RP)
         case(CTD_IGNORE_MODE)
            call multiphase % SetStarMobility(multiphase % M0)
         end select

         call ComputeNSTimeDerivative(mesh, time)

         call multiphase % SetStarMobility(multiphase % M0)

         end select
!
!        -------------------------------------------------------------------------------
!        If IMEX_Explicit, compute cDot with the explicit part of the chemical potential
!        -------------------------------------------------------------------------------
!
         select case (mode)
         case(CTD_IMEX_EXPLICIT)
!$omp do schedule(runtime)
            do eID = 1, size(mesh % elements)
!
!            + Linear part
               mesh % elements(eID) % storage % mu = - IMEX_S0 * mesh % elements(eID) % storage % c &
                                                     - 1.5_RP*(1.0_RP - IMEX_K0)*multiphase % sigma*multiphase % eps*mesh % elements(eID) % storage % cDot
               !mesh % elements(eID) % storage % mu = - IMEX_S0 * mesh % elements(eID) % storage % c &
               !                                      - (1.0_RP - IMEX_K0)* POW2(multiphase % eps)*mesh % elements(eID) % storage % cDot
!
!            + NonLinear part
               !call AddQuarticDWPDerivative(mesh % elements(eID) % storage % c, mesh % elements(eID) % storage % mu)
               call Multiphase_AddChemFEDerivative(mesh % elements(eID) % storage % c, mesh % elements(eID) % storage % mu)
            end do
!$omp end do         
!
!           -----------------------------------
!           Prolong chemical potential to faces
!           -----------------------------------
!
!$omp single
            call mesh % SetStorageToEqn(MU_BC)
            call SetBoundaryConditionsEqn(MU_BC)
!$omp end single
            call mesh % ProlongSolutionToFaces(NCOMP)
!
!           ------------------------------------------------------------
!           Get concentration (lifted) gradients (also prolong to faces)
!           ------------------------------------------------------------
!
            call CHDiscretization % ComputeGradient(NCOMP, NCOMP, mesh, time, chGradientVariables)
!
!           --------------------------------
!           Get chemical potential laplacian
!           --------------------------------
!
!           Get the concentration Laplacian (into QDot => cDot)

            call ComputeLaplacian(mesh, time)

!$omp single
            call mesh % SetStorageToEqn(NS_BC)
            call SetBoundaryConditionsEqn(NS_BC)
!$omp end single
!
!           -----------------------------------------
!           Add the Chemical potential to the NS QDot
!           -----------------------------------------
!
!$omp do schedule(runtime)
            do eID = 1, size(mesh % elements)
               mesh % elements(eID) % storage % QDot(IMC,:,:,:) =   mesh % elements(eID) % storage % QDot(IMC,:,:,:) &
                                                                  + multiphase % M0*mesh % elements(eID) % storage % cDot(1,:,:,:)
            end do
!$omp end do
         end select
!$omp end parallel
!
      END SUBROUTINE ComputeTimeDerivative
!
!////////////////////////////////////////////////////////////////////////////////////
!
!           Navier--Stokes procedures
!           -------------------------
!
!////////////////////////////////////////////////////////////////////////////////////
!
      subroutine ComputeNSTimeDerivative( mesh , t )
         implicit none
         type(HexMesh)              :: mesh
         real(kind=RP)              :: t
         procedure(UserDefinedSourceTermNS_f) :: UserDefinedSourceTermNS
!
!        ---------------
!        Local variables
!        ---------------
!
         integer     :: eID , i, j, k, ierr, fID
         real(kind=RP) :: sqrtRho, invSqrtRho
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
               CALL computeElementInterfaceFlux_MU( f ) 
 
            case (HMESH_BOUNDARY) 
               CALL computeBoundaryFlux_MU(f, t) 
 
            end select 
            end associate 
         end do 
!$omp end do 
!
!        *************************************************************************************
!        Element without shared faces: Surface integrals, scaling of elements with Jacobian, 
!                                      sqrt(rho), and add source terms
!        *************************************************************************************
! 
!$omp do schedule(runtime) private(i,j,k,sqrtRho,invSqrtRho)
         do eID = 1, size(mesh % elements) 
            associate(e => mesh % elements(eID)) 
            if ( e % hasSharedFaces ) cycle
            call TimeDerivative_FacesContribution(e, t, mesh) 
 
            do k = 0, e % Nxyz(3) ; do j = 0, e % Nxyz(2) ; do i = 0, e % Nxyz(1) 
               sqrtRho = sqrt(e % storage % rho(i,j,k))
               invSqrtRho = 1.0_RP / sqrtRho
!
!            + Scale with Jacobian and sqrt(Rho)
               e % storage % QDot(:,i,j,k) = e % storage % QDot(:,i,j,k) * e % geom % InvJacobian(i,j,k)
               e % storage % QDot(IMSQRHOU:IMSQRHOW,i,j,k) = e % storage % QDot(IMSQRHOU:IMSQRHOW,i,j,k) * invSqrtRho
!
!            + Add gravity
               e % storage % QDot(IMSQRHOU:IMSQRHOW,i,j,k) =   e % storage % QDot(IMSQRHOU:IMSQRHOW,i,j,k) & 
                                                             + sqrtRho * dimensionless % invFr2 * dimensionless % gravity_dir

!
!            + Add user defined source terms
               call UserDefinedSourceTermNS(e % geom % x(:,i,j,k), e % storage % Q(:,i,j,k), t, e % storage % S_NS(:,i,j,k), thermodynamics, dimensionless, refValues, multiphase)
               e % storage % QDot(:,i,j,k) = e % storage % QDot(:,i,j,k) + e % storage % S_NS(:,i,j,k) / [1.0_RP,sqrtRho,sqrtRho,sqrtRho,1.0_RP]

            end do         ; end do          ; end do 

            end associate 
         end do
!$omp end do
!
!        ******************************************
!        Do the same for elements with shared faces
!        ******************************************
!
#ifdef _HAS_MPI_
         if ( MPI_Process % doMPIAction ) then
!$omp single
            call mesh % GatherMPIFacesGradients(NCONS)
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
                  CALL computeMPIFaceFlux_MU( f )
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
                  sqrtRho = sqrt(e % storage % rho(i,j,k))
                  invSqrtRho = 1.0_RP / sqrtRho
!   
!               + Scale with Jacobian and sqrt(Rho)
                  e % storage % QDot(:,i,j,k) = e % storage % QDot(:,i,j,k) * e % geom % InvJacobian(i,j,k)
                  e % storage % QDot(IMSQRHOU:IMSQRHOW,i,j,k) = e % storage % QDot(IMSQRHOU:IMSQRHOW,i,j,k) * invSqrtRho
!   
!               + Add gravity
                  e % storage % QDot(IMSQRHOU:IMSQRHOW,i,j,k) =   e % storage % QDot(IMSQRHOU:IMSQRHOW,i,j,k) & 
                                                                + sqrtRho * dimensionless % invFr2 * dimensionless % gravity_dir
   
!   
!               + Add user defined source terms
                  call UserDefinedSourceTermNS(e % geom % x(:,i,j,k), e % storage % Q(:,i,j,k), t, e % storage % S_NS(:,i,j,k), thermodynamics, dimensionless, refValues, multiphase)
                  e % storage % QDot(:,i,j,k) = e % storage % QDot(:,i,j,k) + e % storage % S_NS(:,i,j,k) / [1.0_RP,sqrtRho,sqrtRho,sqrtRho,1.0_RP]
   
               end do         ; end do          ; end do 


               do k = 0, e % Nxyz(3) ; do j = 0, e % Nxyz(2) ; do i = 0, e % Nxyz(1)
                  e % storage % QDot(:,i,j,k) = e % storage % QDot(:,i,j,k) / e % geom % jacobian(i,j,k)
               end do         ; end do          ; end do
               end associate
            end do
!$omp end do
!
!           Add an MPI Barrier
!           ------------------
!$omp single
            call mpi_barrier(MPI_COMM_WORLD, ierr)
!$omp end single
         end if
#endif


      end subroutine ComputeNSTimeDerivative
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
         call HyperbolicDiscretization % ComputeInnerFluxes ( e , mEulerFlux, inviscidContravariantFlux ) 
!
!        Compute viscous contravariant flux
!        ----------------------------------
         call ViscousDiscretization  % ComputeInnerFluxes ( NCONS, NCONS, mViscousFlux, GetmTwoFluidsViscosity, e , viscousContravariantFlux) 
!
!        ************************
!        Perform volume integrals
!        ************************
!
!
!        Compute the total Navier-Stokes flux
!        ------------------------------------
         contravariantFlux = inviscidContravariantFlux - viscousContravariantFlux 
!
!        Perform the Weak Volume Green integral
!        --------------------------------------
         e % storage % QDot = e % storage % QDot + ScalarWeakIntegrals % StdVolumeGreen ( e, NCONS, contravariantFlux ) 


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
      SUBROUTINE computeElementInterfaceFlux_MU(f)
         use FaceClass
         use RiemannSolvers_MU
         IMPLICIT NONE
         TYPE(Face)   , INTENT(inout) :: f   
         integer       :: i, j
         real(kind=RP) :: inv_fluxL(1:NCONS,0:f % Nf(1),0:f % Nf(2))
         real(kind=RP) :: inv_fluxR(1:NCONS,0:f % Nf(1),0:f % Nf(2))
         real(kind=RP) :: visc_flux(1:NCONS,0:f % Nf(1),0:f % Nf(2))
         real(kind=RP) :: fluxL(1:NCONS,0:f % Nf(1),0:f % Nf(2))
         real(kind=RP) :: fluxR(1:NCONS,0:f % Nf(1),0:f % Nf(2))
         real(kind=RP) :: muL, muR
         real(kind=RP) :: UxL(1:NGRAD), UyL(1:NGRAD), UzL(1:NGRAD)
         real(kind=RP) :: UxR(1:NGRAD), UyR(1:NGRAD), UzR(1:NGRAD)

         DO j = 0, f % Nf(2)
            DO i = 0, f % Nf(1)

               call GetmTwoFluidsViscosity(f % storage(1) % Q(IMC,i,j), muL)
               call GetmTwoFluidsViscosity(f % storage(2) % Q(IMC,i,j), muR)

!
!            - Premultiply velocity gradients by the viscosity
!              -----------------------------------------------
               UxL = [1.0_RP,muL,muL,muL,1.0_RP]*f % storage(1) % U_x(:,i,j) 
               UyL = [1.0_RP,muL,muL,muL,1.0_RP]*f % storage(1) % U_y(:,i,j) 
               UzL = [1.0_RP,muL,muL,muL,1.0_RP]*f % storage(1) % U_z(:,i,j) 

               UxR = [1.0_RP,muR,muR,muR,1.0_RP]*f % storage(2) % U_x(:,i,j) 
               UyR = [1.0_RP,muR,muR,muR,1.0_RP]*f % storage(2) % U_y(:,i,j) 
               UzR = [1.0_RP,muR,muR,muR,1.0_RP]*f % storage(2) % U_z(:,i,j) 
!      
!              --------------
!              Viscous fluxes
!              --------------
!      
               CALL ViscousDiscretization % RiemannSolver(nEqn = NCONS, nGradEqn = NCONS, &
                                                  EllipticFlux = mViscousFlux, &
                                                  f = f, &
                                                  QLeft = f % storage(1) % Q(:,i,j), &
                                                  QRight = f % storage(2) % Q(:,i,j), &
                                                  U_xLeft = UxL, &
                                                  U_yLeft = UyL, &
                                                  U_zLeft = UzL, &
                                                  U_xRight = UxR, &
                                                  U_yRight = UyR, &
                                                  U_zRight = UzR, &
                                                  mu_left  = [1.0_RP, multiphase % M0_star, 0.0_RP], &
                                                  mu_right = [1.0_RP, multiphase % M0_star, 0.0_RP], &
                                                  nHat = f % geom % normal(:,i,j) , &
                                                  dWall = f % geom % dWall(i,j), &
                                                  sigma = [multiphase % M0_star, 0.0_RP, 0.0_RP, 0.0_RP, 0.0_RP], &
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
                                  rhoL   = f % storage(1) % rho(i,j), &
                                  rhoR   = f % storage(2) % rho(i,j), &
                                  muL    = f % storage(1) % mu(1,i,j), &
                                  muR    = f % storage(2) % mu(1,i,j), &
                                  nHat   = f % geom % normal(:,i,j), &
                                  t1     = f % geom % t1(:,i,j), &
                                  t2     = f % geom % t2(:,i,j), &
                                  fL     = inv_fluxL(:,i,j), &
                                  fR     = inv_fluxR(:,i,j) )

               
!
!              Multiply by the Jacobian
!              ------------------------
               fluxL(:,i,j) = ( inv_fluxL(:,i,j) - visc_flux(:,i,j)) * f % geom % jacobian(i,j)
               fluxR(:,i,j) = ( inv_fluxR(:,i,j) - visc_flux(:,i,j)) * f % geom % jacobian(i,j)
               
            END DO   
         END DO  
!
!        ---------------------------
!        Return the flux to elements
!        ---------------------------
!
         call f % ProjectFluxToElements(NCONS, fluxL, (/1, HMESH_NONE/))
         call f % ProjectFluxToElements(NCONS, fluxR, (/2, HMESH_NONE/))

      END SUBROUTINE computeElementInterfaceFlux_MU

      SUBROUTINE computeMPIFaceFlux_MU(f)
         use FaceClass
         use RiemannSolvers_MU
         TYPE(Face)   , INTENT(inout) :: f   
         integer       :: i, j
         integer       :: thisSide
         real(kind=RP) :: inv_fluxL(1:NCONS,0:f % Nf(1),0:f % Nf(2))
         real(kind=RP) :: inv_fluxR(1:NCONS,0:f % Nf(1),0:f % Nf(2))
         real(kind=RP) :: visc_flux(1:NCONS,0:f % Nf(1),0:f % Nf(2))
         real(kind=RP) :: fluxL(1:NCONS,0:f % Nf(1),0:f % Nf(2))
         real(kind=RP) :: fluxR(1:NCONS,0:f % Nf(1),0:f % Nf(2))
         real(kind=RP) :: flux(1:NCONS,0:f % Nf(1),0:f % Nf(2),2)
         real(kind=RP) :: muL, muR
         real(kind=RP) :: UxL(1:NGRAD), UyL(1:NGRAD), UzL(1:NGRAD)
         real(kind=RP) :: UxR(1:NGRAD), UyR(1:NGRAD), UzR(1:NGRAD)

         DO j = 0, f % Nf(2)
            DO i = 0, f % Nf(1)

               call GetmTwoFluidsViscosity(f % storage(1) % Q(IMC,i,j), muL)
               call GetmTwoFluidsViscosity(f % storage(2) % Q(IMC,i,j), muR)

!
!            - Premultiply velocity gradients by the viscosity
!              -----------------------------------------------
               UxL = [1.0_RP,muL,muL,muL,1.0_RP]*f % storage(1) % U_x(:,i,j) 
               UyL = [1.0_RP,muL,muL,muL,1.0_RP]*f % storage(1) % U_y(:,i,j) 
               UzL = [1.0_RP,muL,muL,muL,1.0_RP]*f % storage(1) % U_z(:,i,j) 

               UxR = [1.0_RP,muR,muR,muR,1.0_RP]*f % storage(2) % U_x(:,i,j) 
               UyR = [1.0_RP,muR,muR,muR,1.0_RP]*f % storage(2) % U_y(:,i,j) 
               UzR = [1.0_RP,muR,muR,muR,1.0_RP]*f % storage(2) % U_z(:,i,j) 
!      
!              --------------
!              Viscous fluxes
!              --------------
!      
               CALL ViscousDiscretization % RiemannSolver(nEqn = NCONS, nGradEqn = NCONS, &
                                                  EllipticFlux = mViscousFlux, &
                                                  f = f, &
                                                  QLeft = f % storage(1) % Q(:,i,j), &
                                                  QRight = f % storage(2) % Q(:,i,j), &
                                                  U_xLeft = UxL, &
                                                  U_yLeft = UyL, &
                                                  U_zLeft = UzL, &
                                                  U_xRight = UxR, &
                                                  U_yRight = UyR, &
                                                  U_zRight = UzR, &
                                                  mu_left  = [1.0_RP, multiphase % M0_star, 0.0_RP], &
                                                  mu_right = [1.0_RP, multiphase % M0_star, 0.0_RP], &
                                                  nHat = f % geom % normal(:,i,j) , &
                                                  dWall = f % geom % dWall(i,j), &
                                                  sigma = [multiphase % M0_star, 0.0_RP, 0.0_RP, 0.0_RP, 0.0_RP], &
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
                                  rhoL   = f % storage(1) % rho(i,j), &
                                  rhoR   = f % storage(2) % rho(i,j), &
                                  muL    = f % storage(1) % mu(1,i,j), &
                                  muR    = f % storage(2) % mu(1,i,j), &
                                  nHat   = f % geom % normal(:,i,j), &
                                  t1     = f % geom % t1(:,i,j), &
                                  t2     = f % geom % t2(:,i,j), &
                                  fL     = inv_fluxL(:,i,j), &
                                  fR     = inv_fluxR(:,i,j) )

               
!
!              Multiply by the Jacobian
!              ------------------------
               fluxL(:,i,j) = ( inv_fluxL(:,i,j) - visc_flux(:,i,j)) * f % geom % jacobian(i,j)
               fluxR(:,i,j) = ( inv_fluxR(:,i,j) - visc_flux(:,i,j)) * f % geom % jacobian(i,j)
               
            END DO   
         END DO  
!
!        ---------------------------
!        Return the flux to elements
!        ---------------------------
!
         thisSide = maxloc(f % elementIDs, dim = 1)
         flux(:,:,:,1) = fluxL
         flux(:,:,:,2) = fluxR

         call f % ProjectFluxToElements(NCONS, flux(:,:,:,thisSide), (/thisSide, HMESH_NONE/))

      end subroutine ComputeMPIFaceFlux_MU

      SUBROUTINE computeBoundaryFlux_MU(f, time)
      USE ElementClass
      use FaceClass
      USE RiemannSolvers_MU
      IMPLICIT NONE
!
!     ---------
!     Arguments
!     ---------
!
      type(Face),    intent(inout) :: f
      REAL(KIND=RP)                :: time
!
!     ---------------
!     Local variables
!     ---------------
!
      INTEGER                         :: i, j
      INTEGER, DIMENSION(2)           :: N
      REAL(KIND=RP)                   :: inv_fluxL(NCONS), inv_fluxR(NCONS), fv_3d(NCONS,NDIM)
      real(kind=RP)                   :: visc_flux(NCONS, 0:f % Nf(1), 0:f % Nf(2))
      real(kind=RP)                   :: fStar(NCONS, 0:f % Nf(1), 0: f % Nf(2))
      real(kind=RP)                   :: mu
!
!     -------------------
!     Get external states
!     -------------------
!
      do j = 0, f % Nf(2)  ; do i = 0, f % Nf(1)
         f % storage(2) % Q(:,i,j) = f % storage(1) % Q(:,i,j)
         f % storage(2) % mu(1,i,j) = f % storage(1) % mu(1,i,j)
         CALL BCs(f % zone) % bc % FlowState( &
                                      f % geom % x(:,i,j), &
                                      time, &
                                      f % geom % normal(:,i,j), &
                                      f % storage(2) % Q(:,i,j))

         f % storage(2) % rho(i,j) = dimensionless % rho(2) + (dimensionless % rho(1)-dimensionless % rho(2))*f % storage(2) % Q(IMC,i,j)
         f % storage(2) % rho(i,j) = min(max(f % storage(2) % rho(i,j), dimensionless % rho_min),dimensionless % rho_max)
!   
!        --------------
!        Viscous fluxes
!        --------------
!   
         call GetmTwoFluidsViscosity(f % storage(1) % Q(IMC,i,j), mu)

         call mViscousFlux(NCONS, NCONS, f % storage(1) % Q(:,i,j), &
                                         f % storage(1) % U_x(:,i,j), &
                                         f % storage(1) % U_y(:,i,j), &
                                         f % storage(1) % U_z(:,i,j), &
                                         mu, multiphase % M0_star, 0.0_RP, fv_3d)

         visc_flux(:,i,j) =   fv_3d(:,IX)*f % geom % normal(IX,i,j) &
                            + fv_3d(:,IY)*f % geom % normal(IY,i,j) &
                            + fv_3d(:,IZ)*f % geom % normal(IZ,i,j) 


         CALL BCs(f % zone) % bc % FlowNeumann(&
                                           f % geom % x(:,i,j), &
                                           time, &
                                           f % geom % normal(:,i,j), &
                                           f % storage(1) % Q(:,i,j), &
                                           f % storage(1) % U_x(:,i,j), &
                                           f % storage(1) % U_y(:,i,j), &
                                           f % storage(1) % U_z(:,i,j), visc_flux(:,i,j))

!
!           Hyperbolic part
!           -------------
            CALL RiemannSolver(QLeft  = f % storage(1) % Q(:,i,j), &
                               QRight = f % storage(2) % Q(:,i,j), &
                               rhoL   = f % storage(1) % rho(i,j), &
                               rhoR   = f % storage(2) % rho(i,j), &
                               muL    = f % storage(1) % mu(1,i,j), &
                               muR    = f % storage(2) % mu(1,i,j), &
                               nHat   = f % geom % normal(:,i,j), &
                               t1     = f % geom % t1(:,i,j), &
                               t2     = f % geom % t2(:,i,j), &
                               fL     = inv_fluxL, &
                               fR     = inv_fluxR )

            fStar(:,i,j) = (inv_fluxL - visc_flux(:,i,j) ) * f % geom % jacobian(i,j)
         end do   
      end do   

      call f % ProjectFluxToElements(NCONS, fStar, (/1, HMESH_NONE/))

      END SUBROUTINE computeBoundaryFlux_MU
!
!////////////////////////////////////////////////////////////////////////////////////////
!
!        Laplacian procedures
!        --------------------
!
!////////////////////////////////////////////////////////////////////////////////////////
!
      subroutine ComputeLaplacian( mesh , t)
         implicit none
         type(HexMesh)              :: mesh
         real(kind=RP)              :: t
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
               CALL Laplacian_computeElementInterfaceFlux( f ) 
 
            case (HMESH_BOUNDARY) 
               CALL Laplacian_computeBoundaryFlux(f, t) 
            
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
!        ***********************************************************
!        Surface integrals and scaling of elements with shared faces
!        ***********************************************************
! 
#ifdef _HAS_MPI_
         if ( MPI_Process % doMPIAction ) then
!$omp single
            call mesh % GatherMPIFacesGradients(NCOMP)
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
                  CALL Laplacian_computeMPIFaceFlux( f )
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

      end subroutine ComputeLaplacian

      subroutine ComputeLaplacianNeumannBCs( mesh , t)
         implicit none
         type(HexMesh)              :: mesh
         real(kind=RP)              :: t
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
               CALL Laplacian_computeBoundaryFlux(f, t) 
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
         call CHDiscretization  % ComputeInnerFluxes (NCOMP, NCOMP, CHDivergenceFlux, GetCHViscosity, e , contravariantFlux  ) 
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
      subroutine Laplacian_computeElementInterfaceFlux(f)
         use FaceClass
         use Physics
         use PhysicsStorage
         IMPLICIT NONE
         TYPE(Face)   , INTENT(inout) :: f   
         integer       :: i, j
         real(kind=RP) :: flux(1:NCOMP,0:f % Nf(1),0:f % Nf(2))
         real(kind=RP) :: mu

         DO j = 0, f % Nf(2)
            DO i = 0, f % Nf(1)

!      
!              --------------
!              Viscous fluxes
!              --------------
!      
               call GetCHViscosity(0.0_RP, mu)
               CALL CHDiscretization % RiemannSolver(nEqn = NCOMP, nGradEqn = NCOMP, &
                                                  EllipticFlux = CHDivergenceFlux, &
                                                  f = f, &
                                                  QLeft = f % storage(1) % Q(:,i,j), &
                                                  QRight = f % storage(2) % Q(:,i,j), &
                                                  U_xLeft = f % storage(1) % U_x(:,i,j), &
                                                  U_yLeft = f % storage(1) % U_y(:,i,j), &
                                                  U_zLeft = f % storage(1) % U_z(:,i,j), &
                                                  U_xRight = f % storage(2) % U_x(:,i,j), &
                                                  U_yRight = f % storage(2) % U_y(:,i,j), &
                                                  U_zRight = f % storage(2) % U_z(:,i,j), &
                                                  mu_left  = [mu, 0.0_RP, 0.0_RP], &
                                                  mu_right = [mu, 0.0_RP, 0.0_RP], &
                                                  nHat = f % geom % normal(:,i,j) , &
                                                  dWall = f % geom % dWall(i,j), &
                                                  sigma = [1.0_RP], &
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

      end subroutine Laplacian_computeElementInterfaceFlux

      subroutine Laplacian_computeMPIFaceFlux(f)
         use FaceClass
         use Physics
         use PhysicsStorage
         IMPLICIT NONE
         TYPE(Face)   , INTENT(inout) :: f   
         integer       :: i, j
         integer       :: thisSide
         real(kind=RP) :: flux(1:NCOMP,0:f % Nf(1),0:f % Nf(2))
         real(kind=RP) :: mu

         DO j = 0, f % Nf(2)
            DO i = 0, f % Nf(1)
!      
!              --------------
!              Viscous fluxes
!              --------------
!      
               call GetCHViscosity(0.0_RP, mu)
               CALL CHDiscretization % RiemannSolver(nEqn = NCOMP, nGradEqn = NCOMP, &
                                                  EllipticFlux = CHDivergenceFlux, &
                                                  f = f, &
                                                  QLeft = f % storage(1) % Q(:,i,j), &
                                                  QRight = f % storage(2) % Q(:,i,j), &
                                                  U_xLeft = f % storage(1) % U_x(:,i,j), &
                                                  U_yLeft = f % storage(1) % U_y(:,i,j), &
                                                  U_zLeft = f % storage(1) % U_z(:,i,j), &
                                                  U_xRight = f % storage(2) % U_x(:,i,j), &
                                                  U_yRight = f % storage(2) % U_y(:,i,j), &
                                                  U_zRight = f % storage(2) % U_z(:,i,j), &
                                                  mu_left  = [mu, 0.0_RP, 0.0_RP], &
                                                  mu_right = [mu, 0.0_RP, 0.0_RP], &
                                                  nHat = f % geom % normal(:,i,j) , &
                                                  dWall = f % geom % dWall(i,j), &
                                                  sigma = [1.0_RP], &
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

      end subroutine Laplacian_ComputeMPIFaceFlux

      subroutine Laplacian_computeBoundaryFlux(f, time)
      USE ElementClass
      use FaceClass
      USE EllipticDiscretizations
      USE Physics
      use PhysicsStorage
      IMPLICIT NONE
!
!     ---------
!     Arguments
!     ---------
!
      type(Face),    intent(inout) :: f
      REAL(KIND=RP)                :: time
!
!     ---------------
!     Local variables
!     ---------------
!
      INTEGER                         :: i, j
      INTEGER, DIMENSION(2)           :: N
      real(kind=RP)                   :: flux(NCOMP, 0:f % Nf(1), 0:f % Nf(2)), fv_3d(NCOMP,NDIM)
      real(kind=RP)                   :: mu
!
!     -------------------
!     Get external states
!     -------------------
!
      do j = 0, f % Nf(2)  ; do i = 0, f % Nf(1)
!
!        --------------
!        Viscous fluxes
!        --------------
!   
         f % storage(2) % Q(:,i,j) = f % storage(1) % Q(:,i,j)
         call GetCHViscosity(0.0_RP, mu)
         call CHDivergenceFlux(NCOMP, NCOMP, f % storage(1) % Q(:,i,j), &
                                             f % storage(1) % U_x(:,i,j), &
                                             f % storage(1) % U_y(:,i,j), &
                                             f % storage(1) % U_z(:,i,j), &
                                             mu, 0.0_RP, 0.0_RP, fv_3d)

         flux(:,i,j) =   fv_3d(:,IX)*f % geom % normal(IX,i,j) &
                       + fv_3d(:,IY)*f % geom % normal(IY,i,j) &
                       + fv_3d(:,IZ)*f % geom % normal(IZ,i,j) 

         CALL BCs(f % zone) % bc % NeumannForEqn(NCOMP, NCOMP, &
                                           f % geom % x(:,i,j), &
                                           time, &
                                           f % geom % normal(:,i,j), &
                                           f % storage(1) % Q(:,i,j), &
                                           f % storage(1) % U_x(:,i,j), &
                                           f % storage(1) % U_y(:,i,j), &
                                           f % storage(1) % U_z(:,i,j), flux(:,i,j))

         flux(:,i,j) = flux(:,i,j) * f % geom % jacobian(i,j)

      end do               ; end do

      call f % ProjectFluxToElements(NCOMP, flux, (/1, HMESH_NONE/))

      end subroutine Laplacian_computeBoundaryFlux

      SUBROUTINE ComputeTimeDerivativeIsolated( mesh, particles, time, mode)
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
         integer,             intent(in)  :: mode

      end subroutine ComputeTimeDerivativeIsolated

end module SpatialDiscretization
