!
!//////////////////////////////////////////////////////
!
!   @File:    SpatialDiscretization.f90
!   @Author:  Juan (juan.manzanero@upm.es)
!   @Created: Tue Apr 24 17:10:06 2018
!   @Last revision date: Sun Sep 30 21:41:41 2018
!   @Last revision author: Andrés Rueda (am.rueda@upm.es)
!   @Last revision commit: 6ccda27143afdf4445c53d1d8364e5cff10baabc
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
      use GradientsStabilization
      use BoundaryConditions, only: BCs, SetBoundaryConditionsEqn, NS_BC, C_BC, MU_BC
#ifdef _HAS_MPI_
      use mpi
#endif

      private
      public  ComputeLaplacian, DGSpatial_ComputeGradient
      public  Initialize_SpaceAndTimeMethods, ComputeTimeDerivative, ComputeTimeDerivativeIsolated
      public  ComputeTimeDerivative_onlyLinear, ComputetimeDerivative_onlyNonLinear
      public  Finalize_SpaceAndTimeMethods
      public  viscousDiscretizationKey, CHDiscretizationKey

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

            call ViscousDiscretization % Construct(controlVariables, iViscousFlux0D, iViscousFlux2D, iViscousFlux3D, GetiNSCHViscosity, "NS")
            call ViscousDiscretization % Describe
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
   
            call CHDiscretization % Construct(controlVariables, CHDivergenceFlux0D, CHDivergenceFlux2D, CHDivergenceFlux3D, GetCHViscosity, "CH")
            call CHDiscretization % Describe
         
         end if
!
!        Compute wall distances
!        ----------------------
         call mesh % ComputeWallDistances

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
         integer,             intent(in) :: mode
!
!        ---------------
!        Local variables
!        ---------------
!
         class(Element), pointer    :: e
         INTEGER :: k, eID, fID, i, j
         logical  :: CH_enable_linear, CH_enable_nonlinear, NS_enable, CH_enable
!
!        Configure the time derivative
!        -----------------------------
         select case(mode)
         case (CTD_IGNORE_MODE)
            NS_enable = .true.  ; CH_enable_linear = .true.  ; CH_enable_nonlinear = .true.

         case (CTD_ONLY_NS)
            NS_enable = .true.  ; CH_enable_linear = .false. ; CH_enable_nonlinear = .false.

         case (CTD_NS_AND_CH)
            NS_enable = .true.  ; CH_enable_linear = .true.  ; CH_enable_nonlinear = .true.

         case (CTD_ONLY_CH)
            NS_enable = .false. ; CH_enable_linear = .true.  ; CH_enable_nonlinear = .true.

         case (CTD_ONLY_CH_LIN)
            NS_enable = .false. ; CH_enable_linear = .true.  ; CH_enable_nonlinear = .false.

         case (CTD_ONLY_CH_NONLIN)
            NS_enable = .false. ; CH_enable_linear = .false. ; CH_enable_nonlinear = .true.

         case default
            print*, "Unrecognized mode"
            errorMessage(STD_OUT)
            stop

         end select

         if ( (.not. CH_enable_linear) .and. (.not. CH_enable_nonlinear) ) then
            CH_enable = .false.
         else
            CH_enable = .true.
         end if
         
!
!        -----------------------------------------
!        Prolongation of the solution to the faces
!        -----------------------------------------
!
!$omp parallel shared(mesh, time) private(k, eID, fID, i, j)
!
!        *****************************
!        Obtain the NS time derivative
!        *****************************
!
         if ( NS_enable ) then
!$omp single
            call mesh % SetStorageToEqn(NS_BC)
            call SetBoundaryConditionsEqn(NS_BC)
!$omp end single

            call mesh % ProlongSolutionToFaces(NINC)
!
!           ----------------
!           Update MPI Faces
!           ----------------
!
#ifdef _HAS_MPI_
!$omp single
errorMessage(STD_OUT)
stop
         !call mesh % UpdateMPIFacesSolution
!$omp end single
#endif
!
!           -----------------
!           Compute gradients
!           -----------------
!
            if ( computeGradients ) then
               call ViscousDiscretization % ComputeGradient( NINC, NINC, mesh , time , iNSGradientValuesForQ_0D, iNSGradientValuesForQ_3D)
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
!           -----------------------
!           Compute time derivative
!           -----------------------
!
            call ComputeNSTimeDerivative(mesh              = mesh , &
                                         particles         = particles, &
                                         t                 = time)
         end if
!
!        *****************************************
!        Compute the Cahn-Hilliard time derivative
!        *****************************************
!
!        ------------------------------
!        Change memory to concentration
!        ------------------------------

         if ( CH_enable ) then
!$omp single
            call mesh % SetStorageToEqn(C_BC)
            call SetBoundaryConditionsEqn(C_BC)
!$omp end single
!
!           Prolong solution to faces
!           -------------------------
            call mesh % ProlongSolutionToFaces(NCOMP)
!   
!           ----------------
!           Update MPI Faces
!           ----------------
!   
#ifdef _HAS_MPI_
!$omp single
   errorMessage(STD_OUT)
   stop
            !call mesh % UpdateMPIFacesSolution
!$omp end single
#endif
!   
!           -----------------
!           Compute gradients: prolongation has already been performe: prolongation has already been performedd
!           -----------------
!   
            call CHDiscretization % ComputeGradient( NCOMP, NCOMP, mesh , time, CHGradientValuesForQ_0D, CHGradientValuesForQ_3D)
   
#ifdef _HAS_MPI_
!$omp single
   errorMessage(STD_OUT)
   stop
            !call mesh % UpdateMPIFacesGradients
!$omp end single
#endif
!   
!           ------------------------------
!           Compute the chemical potential
!           ------------------------------
!   
!           Linear part
!           -----------
            if ( CH_enable_linear ) then
               call ComputeLaplacian(mesh = mesh , &
                                     t    = time)
            else
       
               call ComputeLaplacianNeumannBCs(mesh = mesh , &
                                     t    = time)
            end if

            if ( CH_enable_nonlinear) then
!$omp do schedule(runtime)
               do eID = 1, mesh % no_of_elements
                  e => mesh % elements(eID)
                  e % storage % mu = - POW2(multiphase % eps) * e % storage % QDot
                  call AddQuarticDWPDerivative(e % storage % c, e % storage % mu)
!   
!                 Move storage to chemical potential
!                 ----------------------------------
                  call e % storage % SetStorageToCH_mu
               end do
!$omp end do
            else
               do eID = 1, mesh % no_of_elements
                  e => mesh % elements(eID)
                  e % storage % mu = - POW2(multiphase % eps) * e % storage % QDot
!   
!                 Move storage to chemical potential
!                 ----------------------------------
                  call e % storage % SetStorageToCH_mu
               end do
            end if
   
   
!$omp single
            call mesh % SetStorageToEqn(MU_BC)
            call SetBoundaryConditionsEqn(MU_BC)
!$omp end single
!   
!           *************************
!           Compute cDot: Q stores mu
!           *************************
!   
!           -----------------------------------------
!           Prolongation of the solution to the faces
!           -----------------------------------------
!   
            call mesh % ProlongSolutionToFaces(NCOMP)
!   
!           ----------------
!           Update MPI Faces
!           ----------------
!   
#ifdef _HAS_MPI_
!$omp single
   errorMessage(STD_OUT)
   stop
            !call mesh % UpdateMPIFacesSolution
!$omp end single
#endif
!   
!           -----------------
!           Compute gradients
!           -----------------
!   
            call CHDiscretization % ComputeGradient( NCOMP, NCOMP, mesh , time, CHGradientValuesForQ_0D, CHGradientValuesForQ_3D)
   
#ifdef _HAS_MPI_
!$omp single
!            call mesh % UpdateMPIFacesGradients
   errorMessage(STD_OUT)
   stop
!$omp end single
#endif
!   
!           ------------------------------
!           Compute the chemical potential
!           ------------------------------
!   
            call ComputeLaplacian(mesh = mesh , &
                                  t    = time )
!   
!           Scale QDot with the Peclet number
!           ---------------------------------
!$omp do schedule(runtime) private(e)
            do eID = 1, mesh % no_of_elements
               e => mesh % elements(eID)
               e % storage % QDot = (1.0_RP / multiphase % Pe) * e % storage % QDot
            end do
!$omp end do
!   
!           *****************************
!           Return the concentration to Q
!           *****************************
!   
!$omp single
            call mesh % SetStorageToEqn(C_BC)
            call SetBoundaryConditionsEqn(C_BC)
!$omp end single

         end if  ! CH_enable
!
!        ****************************
!        Return NS as default storage
!        ****************************
!
         if ( NS_enable ) then
!$omp single
            call mesh % SetStorageToEqn(NS_BC)
            call SetBoundaryConditionsEqn(NS_BC)
!$omp end single
!
!           ****************************
!           Compute the Capilar pressure
!           ****************************
!
!$omp do schedule(runtime) private(e,i,j,k)
            do eID = 1, mesh % no_of_elements
               e => mesh % elements(eID) 
               do k = 0, e % Nxyz(3) ; do j = 0, e % Nxyz(2)   ; do i = 0, e % Nxyz(1)
                  e % storage % QDot(INSRHOU,i,j,k) =   e % storage % QDot(INSRHOU,i,j,k) &
                                                         + (1.0_RP / (multiphase % eps * dimensionless % Re * multiphase % Ca)) * e % storage % mu(1,i,j,k) * e % storage % c_x(1,i,j,k) 
                  e % storage % QDot(INSRHOV,i,j,k) =   e % storage % QDot(INSRHOV,i,j,k) &
                                                         + (1.0_RP / (multiphase % eps * dimensionless % Re * multiphase % Ca)) * e % storage % mu(1,i,j,k) * e % storage % c_y(1,i,j,k)  
                  e % storage % QDot(INSRHOW,i,j,k) =   e % storage % QDot(INSRHOW,i,j,k) &
                                                         + (1.0_RP / (multiphase % eps * dimensionless % Re * multiphase % Ca)) * e % storage % mu(1,i,j,k) * e % storage % c_z(1,i,j,k)
   
               end do                ; end do                  ; end do
            end do
!$omp end do

         end if ! NS_enable
!$omp end parallel
!
      END SUBROUTINE ComputeTimeDerivative
!
!////////////////////////////////////////////////////////////////////////
!
!     This routine computes the time derivative element by element, without considering the Riemann Solvers
!     This is useful for estimating the isolated truncation error
!
      SUBROUTINE ComputeTimeDerivativeIsolated( mesh, particles, time, mode)
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
         integer,             intent(in) :: mode
      END SUBROUTINE ComputeTimeDerivativeIsolated

      subroutine ComputeNSTimeDerivative( mesh , particles, t)
         implicit none
         type(HexMesh)              :: mesh
         type(Particles_t)          :: particles
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
               CALL computeBoundaryFlux_iNS(f, t) 
 
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
            call TimeDerivative_FacesContribution(e, t, mesh,NINC) 
 
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
               call TimeDerivative_FacesContribution(e, t, mesh, NINC) 
   
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
                  e % storage % QDot(INSRHOU:INSRHOW,i,j,k) = e % storage % QDot(INSRHOU:INSRHOW,i,j,k) + &
                                                        e % storage % Q(INSRHO,i,j,k) * &
                                    dimensionless % invFr2 * dimensionless % gravity_dir

               end do                ; end do                ; end do
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
      subroutine ComputeNSTimeDerivativeIsolated( mesh , t )
         implicit none
         type(HexMesh)              :: mesh
         real(kind=RP)              :: t
      end subroutine ComputeNSTimeDerivativeIsolated
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
      subroutine TimeDerivative_FacesContribution( e , t , mesh, nEqn)
         use HexMeshClass
         implicit none
         type(Element)           :: e
         real(kind=RP)           :: t
         type(HexMesh)           :: mesh
         integer :: nEqn

         e % storage % QDot = e % storage % QDot - ScalarWeakIntegrals % StdFace( e, nEqn, &
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

               call ViscousDiscretization % GetViscosity(f % storage(1) % c(1,i,j), muL)
               call ViscousDiscretization % GetViscosity(f % storage(2) % c(1,i,j), muR)
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
                                                  mu   = mu, beta = 0.0_RP, kappa = 0.0_RP, &
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
               call ViscousDiscretization % GetViscosity(f % storage(1) % c(1,i,j), mu)

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
                                                  mu   = mu, beta = 0.0_RP, kappa = 0.0_RP, &
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

      SUBROUTINE computeBoundaryFlux_iNS(f, time)
      USE ElementClass
      use FaceClass
      USE RiemannSolvers_iNS
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
      REAL(KIND=RP)                   :: inv_flux(NINC)
      real(kind=RP)                   :: visc_flux(NINC, 0:f % Nf(1), 0:f % Nf(2))
      real(kind=RP)                   :: fStar(NINC, 0:f % Nf(1), 0: f % Nf(2))
      real(kind=RP)                   :: mu
      
!
!     -------------------
!     Get external states
!     -------------------
!
      do j = 0, f % Nf(2)  ; do i = 0, f % Nf(1)
         f % storage(2) % Q(:,i,j) = f % storage(1) % Q(:,i,j)
         CALL BCs(f % zone) % bc % FlowState( &
                                      f % geom % x(:,i,j), &
                                      time, &
                                      f % geom % normal(:,i,j), &
                                      f % storage(2) % Q(:,i,j))

      end do               ; end do

      if ( .true. ) then
         do j = 0, f % Nf(2)  ; do i = 0, f % Nf(1)
            f % storage(2) % U_x(:,i,j) = f % storage(1) % U_x(:,i,j)
            f % storage(2) % U_y(:,i,j) = f % storage(1) % U_y(:,i,j)
            f % storage(2) % U_z(:,i,j) = f % storage(1) % U_z(:,i,j)

            CALL BCs(f % zone) % bc % FlowNeumann(&
                                              f % geom % x(:,i,j), &
                                              time, &
                                              f % geom % normal(:,i,j), &
                                              f % storage(2) % Q(:,i,j), &
                                              f % storage(2) % U_x(:,i,j), &
                                              f % storage(2) % U_y(:,i,j), &
                                              f % storage(2) % U_z(:,i,j))

!   
!           --------------
!           Viscous fluxes
!           --------------
!   
            call ViscousDiscretization % GetViscosity(f % storage(1) % c(1,i,j), mu)

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
                                               mu   = mu, beta = 0.0_RP, kappa = 0.0_RP, &
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
!///////////////////////////////////////////////////////////////////////////////////////////
!
!           Procedures to compute the state variables Laplacian
!           ---------------------------------------------------
!
!///////////////////////////////////////////////////////////////////////////////////////////
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
               CALL computeElementInterfaceFlux_Laplacian( f ) 
 
            case (HMESH_BOUNDARY) 
               CALL computeBoundaryFlux_Laplacian(f, t) 
            
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
               CALL computeBoundaryFlux_Laplacian(f, t) 
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
            call TimeDerivative_FacesContribution(e, t, mesh, NCOMP) 
 
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
               call TimeDerivative_FacesContribution(e, t, mesh,NCOMP) 
 
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
                                                  mu   = 1.0_RP, beta = 0.0_RP, kappa = 0.0_RP, &
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
                                                  mu   = 1.0_RP, beta = 0.0_RP, kappa = 0.0_RP, &
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

      subroutine computeBoundaryFlux_Laplacian(f, time)
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
      real(kind=RP)                   :: flux(NCOMP, 0:f % Nf(1), 0:f % Nf(2))
!
!     -------------------
!     Get external states
!     -------------------
!
      do j = 0, f % Nf(2)  ; do i = 0, f % Nf(1)
         f % storage(2) % Q(:,i,j) = f % storage(1) % Q(:,i,j)
         CALL BCs(f % zone) % bc % StateForEqn(NCOMP, f % geom % x(:,i,j), &
                                      time, &
                                      f % geom % normal(:,i,j), &
                                      f % storage(2) % Q(:,i,j))

      end do               ; end do

      do j = 0, f % Nf(2)  ; do i = 0, f % Nf(1)
         f % storage(2) % U_x(:,i,j) = f % storage(1) % U_x(:,i,j)
         f % storage(2) % U_y(:,i,j) = f % storage(1) % U_y(:,i,j)
         f % storage(2) % U_z(:,i,j) = f % storage(1) % U_z(:,i,j)

         CALL BCs(f % zone) % bc % NeumannForEqn(NCOMP, NCOMP, &
                                           f % geom % x(:,i,j), &
                                           time, &
                                           f % geom % normal(:,i,j), &
                                           f % storage(2) % Q(:,i,j), &
                                           f % storage(2) % U_x(:,i,j), &
                                           f % storage(2) % U_y(:,i,j), &
                                           f % storage(2) % U_z(:,i,j))

         f % storage(1) % U_x(:,i,j) = f % storage(2) % U_x(:,i,j)
         f % storage(1) % U_y(:,i,j) = f % storage(2) % U_y(:,i,j)
         f % storage(1) % U_z(:,i,j) = f % storage(2) % U_z(:,i,j)
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
                                            mu = 1.0_RP, beta = 0.0_RP, kappa = 0.0_RP, &
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
