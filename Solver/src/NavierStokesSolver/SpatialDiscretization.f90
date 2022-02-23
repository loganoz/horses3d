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
      use VariableConversion, only: NSGradientVariables_STATE, GetNSViscosity, NSGradientVariables_ENTROPY, &
                                    GetGradientValues_f, NSGradientVariables_ENERGY, get_laminar_mu_kappa, &
                                    set_getVelocityGradients
      use ProblemFileFunctions, only: UserDefinedSourceTermNS_f
      use BoundaryConditions
#ifdef _HAS_MPI_
      use mpi
#endif
	USE OMP_LIB

      private
      public   ComputeTimeDerivative, ComputeTimeDerivativeIsolated, viscousDiscretizationKey
      public   Initialize_SpaceAndTimeMethods, Finalize_SpaceAndTimeMethods

!!!!!!!
! DEFINE TURBINE, BLADE, AIRFOIL
!!!!!!!

    type airfoil
    integer                  :: num_aoa 
    real(KIND=RP), allocatable      :: aoa(:)  ! in rad
    real(KIND=RP), allocatable      :: cl(:)  ! in rad
    real(KIND=RP), allocatable      :: cd(:)  ! in rad
    end type

    type blade
    real(KIND=RP), allocatable      :: r_R(:)  ! in mm
    real(KIND=RP), allocatable      :: chord(:)   ! in mm
    real(KIND=RP), allocatable      :: twist(:)    ! in rad
    real(KIND=RP), allocatable      :: azimuth_angle(:)   ! in rad
    integer, allocatable     :: num_airfoils(:)    ! in rad
    CHARACTER(LEN=30), allocatable :: airfoil_files(:,:)  ! file names for Cl-Cd
    type(airfoil), allocatable     :: airfoil(:)  ! airfoil data AoA-Cl-Cd
    real(KIND=RP), allocatable      :: local_velocity(:)  ! local flow speed at blade section, im m/s
    real(KIND=RP), allocatable      :: local_angle(:)  ! local flow angle at blade section, in rad
    real(KIND=RP), allocatable      :: local_lift(:)  ! blade sectional Lift
    real(KIND=RP), allocatable      :: local_drag(:)  ! blade sectional Lift
    real(KIND=RP), allocatable      :: point_xyz_loc(:,:) ! x,y location of blade points
    real(KIND=RP), allocatable      :: local_torque(:)  ! Nm
    real(KIND=RP), allocatable      :: local_thrust(:)  ! N
    real(KIND=RP), allocatable      :: local_root_bending(:)  ! Nm
    end type

    type turbine
    integer		           :: num_blades=3  ! number of blades -> harcoded 3 blades
    real(KIND=RP)                  :: radius ! turb radius in mm
    real(KIND=RP)                  :: blade_pitch ! turb radius in rad
    real(KIND=RP)                  :: rot_speed ! rad/s
    real(KIND=RP)                  :: hub_cood_x,hub_cood_y,hub_cood_z ! hub height in mm
    real(KIND=RP)                  :: normal_x, normal_y, normal_z ! rotor normal pointing backward
    type(blade)             	   :: blade(3) ! hardcoded 3 blades
    integer                 	   :: num_blade_sections  ! numer of 2D section for Cl-Cd data
    real(KIND=RP)                  :: blade_torque(3)
    real(KIND=RP)                  :: blade_thrust(3)
    real(KIND=RP)                  :: blade_root_bending(3)
    end type
                                   
    type farm
    integer :: num_turbines
    type(turbine), allocatable   :: turbine(:)
    end type


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

         SUBROUTINE computeBoundaryFluxF(f, time, mesh)
            use SMConstants
            use FaceClass,  only: Face
            use HexMeshClass
            IMPLICIT NONE
            type(Face),    intent(inout) :: f
            REAL(KIND=RP)                :: time
            type(HexMesh), intent(in)   :: mesh
         end subroutine computeBoundaryFluxF
      end interface
      
      procedure(computeElementInterfaceFluxF), pointer :: computeElementInterfaceFlux
      procedure(computeMPIFaceFluxF),          pointer :: computeMPIFaceFlux
      procedure(computeBoundaryFluxF),         pointer :: computeBoundaryFlux

      procedure(GetGradientValues_f),           pointer :: GetGradients
      procedure(EllipticFlux_f),                pointer :: ViscousFlux
   
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
         use WallFunctionConnectivity
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
         character(len=*), parameter      :: gradient_variables_key = "gradient variables"
         character(len=LINE_LENGTH)       :: gradient_variables
         
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
               if ( controlVariables % ContainsKey(gradient_variables_key) ) then   
                  gradient_variables = controlVariables % StringValueForKey(gradient_variables_key, LINE_LENGTH)   
                  call toLower(gradient_variables)

                  select case (trim(gradient_variables))
                  case ("state")
                     call SetGradientVariables(GRADVARS_STATE)
                     GetGradients => NSGradientVariables_STATE
                     ViscousFlux  => ViscousFlux_STATE
                     call set_getVelocityGradients(GRADVARS_STATE)
   
                  case ("entropy")
                     call SetGradientVariables(GRADVARS_ENTROPY)
                     GetGradients => NSGradientVariables_ENTROPY
                     ViscousFlux  => ViscousFlux_ENTROPY
                     call set_getVelocityGradients(GRADVARS_ENTROPY)

                  case ("energy")
                     call SetGradientVariables(GRADVARS_ENERGY)
                     GetGradients => NSGradientVariables_ENERGY
                     ViscousFlux  => ViscousFlux_ENERGY
                     call set_getVelocityGradients(GRADVARS_ENERGY)

                  case default
                     print*, 'Entropy variables "',trim(gradient_variables),'" are not currently implemented'
                     write(STD_OUT,'(A)') "Implemented options are:"
                     write(STD_OUT,'(A)') "  * State"
                     write(STD_OUT,'(A)') "  * Entropy"
                     write(STD_OUT,'(A)') "  * Energy"
                     errorMessage(STD_OUT)
                     stop
                  end select

               else
!
!                 Set state as default option
!                 ---------------------------
                  call SetGradientVariables(GRADVARS_STATE)
                  GetGradients => NSGradientVariables_STATE
                  ViscousFlux  => ViscousFlux_STATE
                  call set_getVelocityGradients(GRADVARS_STATE)

               end if

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

               call ViscousDiscretization % Construct(controlVariables, ELLIPTIC_NS)
               call ViscousDiscretization % Describe
      
            else
               if (.not. allocated(ViscousDiscretization)) allocate(EllipticDiscretization_t :: ViscousDiscretization)
               call ViscousDiscretization % Construct(controlVariables, ELLIPTIC_NS)
!
!              Set state as default option
!              ---------------------------
               call SetGradientVariables(GRADVARS_STATE)
               GetGradients => NSGradientVariables_STATE
               ViscousFlux  => ViscousFlux_STATE
               call set_getVelocityGradients(GRADVARS_STATE)
               
            end if

   !
   !        Initialize models
   !        -----------------
            call InitializeLESModel(LESModel, controlVariables)
         
         end if
!
!        Initialize SVV
!        --------------
         call InitializeSVV(SVV, controlVariables, mesh)
         
         if (.not. mesh % child) then
            if ( SVV % enabled ) then
               computeElementInterfaceFlux => computeElementInterfaceFlux_SVV
               computeMPIFaceFlux          => computeMPIFaceFlux_SVV
               computeBoundaryFlux         => computeBoundaryFlux_SVV
            else
               computeElementInterfaceFlux => computeElementInterfaceFlux_NS
               computeMPIFaceFlux          => computeMPIFaceFlux_NS
               computeBoundaryFlux         => computeBoundaryFlux_NS
            end if
         end if

!
!        Initialize Wall Function
!        ------------------------
         call Initialize_WallConnection(controlVariables, mesh)
         
      end subroutine Initialize_SpaceAndTimeMethods
!
!////////////////////////////////////////////////////////////////////////
!
      subroutine Finalize_SpaceAndTimeMethods
         implicit none
         IF ( ALLOCATED(HyperbolicDiscretization) ) DEALLOCATE( HyperbolicDiscretization )
         IF ( ALLOCATED(LESModel) )       DEALLOCATE( LESModel )
         
         call SVV % destruct()
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
         INTEGER :: k

         call SetBoundaryConditionsEqn(NS_BC)
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
#ifdef _HAS_MPI_
!$omp single
         call mesh % UpdateMPIFacesSolution(NCONS)
!$omp end single
#endif
!
!        -----------------
!        Compute gradients
!        -----------------
!
         if ( computeGradients ) then
            call ViscousDiscretization % ComputeGradient( NCONS, NGRAD, mesh , time, GetGradients)
         end if

#ifdef _HAS_MPI_
!$omp single
         if ( flowIsNavierStokes ) then
            call mesh % UpdateMPIFacesGradients(NGRAD)
         end if
!$omp end single
#endif
!         call ComputeArtificialViscosity(mesh)
!
!        -----------------------
!        Compute time derivative
!        -----------------------
!
         call TimeDerivative_ComputeQDot(mesh = mesh , &
                                         particles = particles, &
                                         t    = time)
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
         TYPE(HexMesh), target      :: mesh
         type(Particles_t)          :: particles
         REAL(KIND=RP)              :: time
         integer,             intent(in)  :: mode
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
            CALL BaseClass_ComputeGradient( ViscousDiscretization, NCONS, NGRAD, mesh , time , GetGradients)
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

      subroutine TimeDerivative_ComputeQDot( mesh , particles, t)
         implicit none
         type(HexMesh)              :: mesh
         type(Particles_t)          :: particles
         real(kind=RP)              :: t
         procedure(UserDefinedSourceTermNS_f) :: UserDefinedSourceTermNS
!
!        ---------------
!        Local variables
!        ---------------
!
         integer     :: eID , i, j, k, ierr, fID, iFace, iEl
         real(kind=RP)  :: mu_smag, delta

	real(kind=RP) :: density, Cl, Cd, wind_speed, aoa, theta, interp, epsil, Rotor_force, red_factor, &
			t_init,c1,c2, tip_correct,g1_func 

    	integer    ::  io, ii, jj
    	CHARACTER(LEN=40) :: arg, char1

   !!!!!!!!
    type (farm) :: farm1 ! only one farm
!
!        ***********************************************
!        Compute the viscosity at the elements and faces
!        ***********************************************
!
         if (flowIsNavierStokes) then
!$omp do schedule(runtime) private(i,j,k)
            do eID = 1, size(mesh % elements)
               associate(e => mesh % elements(eID))
               do k = 0, e % Nxyz(3) ; do j = 0, e % Nxyz(2) ; do i = 0, e % Nxyz(1)
                  call get_laminar_mu_kappa(e % storage % Q(:,i,j,k), e % storage % mu_NS(1,i,j,k), e % storage % mu_NS(2,i,j,k))
               end do                ; end do                ; end do
               end associate
            end do
!$omp end do
         end if


         if ( LESModel % active) then
!$omp do schedule(runtime) private(i,j,k,delta,mu_smag)
            do eID = 1, size(mesh % elements)
               associate(e => mesh % elements(eID))
               delta = ((e % geom % Volume ) ** (1.0_RP / 3.0_RP))/(real(e % Nxyz (1)) -1.0_RP)
               !delta = (e % geom % Volume / product(e % Nxyz + 1)) ** (1.0_RP / 3.0_RP)
               do k = 0, e % Nxyz(3) ; do j = 0, e % Nxyz(2) ; do i = 0, e % Nxyz(1)
                  call LESModel % ComputeViscosity(delta, e % geom % dWall(i,j,k), e % storage % Q(:,i,j,k),   &
                                                                                   e % storage % U_x(:,i,j,k), &
                                                                                   e % storage % U_y(:,i,j,k), &
                                                                                   e % storage % U_z(:,i,j,k), &
                                                                                   mu_smag)
                  e % storage % mu_NS(1,i,j,k) = e % storage % mu_NS(1,i,j,k) + mu_smag
                  e % storage % mu_NS(2,i,j,k) = e % storage % mu_NS(2,i,j,k) + mu_smag * dimensionless % mu_to_kappa
               end do                ; end do                ; end do
               end associate
            end do
!$omp end do
      end if
!
!        Compute viscosity at interior and boundary faces
!        ------------------------------------------------
         call compute_viscosity_at_faces(size(mesh % faces_interior), 2, mesh % faces_interior, mesh)
         call compute_viscosity_at_faces(size(mesh % faces_boundary), 1, mesh % faces_boundary, mesh)
!
!        ****************
!        Volume integrals
!        ****************
!
!$omp do schedule(runtime) 
         do eID = 1 , size(mesh % elements)
            call TimeDerivative_VolumetricContribution( mesh, mesh % elements(eID) , t)
         end do
!$omp end do 
!
!        ******************************************
!        Compute Riemann solver of non-shared faces
!        ******************************************
!
!$omp do schedule(runtime) private(fID)
         do iFace = 1, size(mesh % faces_interior)
            fID = mesh % faces_interior(iFace)
            call computeElementInterfaceFlux(mesh % faces(fID))
         end do
!$omp end do nowait

!$omp do schedule(runtime) private(fID)
         do iFace = 1, size(mesh % faces_boundary)
            fID = mesh % faces_boundary(iFace)
            call computeBoundaryFlux(mesh % faces(fID), t, mesh)
         end do
!$omp end do 
!
!        ***************************************************************
!        Surface integrals and scaling of elements with non-shared faces
!        ***************************************************************
! 
!$omp do schedule(runtime) private(i,j,k,eID)
         do iEl = 1, size(mesh % elements_sequential)
            eID = mesh % elements_sequential(iEl)
            associate(e => mesh % elements(eID)) 
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
               call mesh % GatherMPIFacesGradients(NGRAD)
            else  
               call mesh % GatherMPIFacesSolution(NCONS)
            end if          
!$omp end single
!
!           Compute viscosity at MPI faces
!           ------------------------------
            call compute_viscosity_at_faces(size(mesh % faces_mpi), 2, mesh % faces_mpi, mesh)
!
!           **************************************
!           Compute Riemann solver of shared faces
!           **************************************
!
!$omp do schedule(runtime) private(fID)
            do iFace = 1, size(mesh % faces_mpi)
               fID = mesh % faces_mpi(iFace)
               call computeMPIFaceFlux(mesh % faces(fID))
            end do
!$omp end do 
!
!           ***********************************************************
!           Surface integrals and scaling of elements with shared faces
!           ***********************************************************
! 
!$omp do schedule(runtime) private(i,j,k,eID)
            do iEl = 1, size(mesh % elements_mpi)
               eID = mesh % elements_mpi(iEl)
               associate(e => mesh % elements(eID)) 
               call TimeDerivative_FacesContribution(e, t, mesh) 
   
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
!


    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!$omp critical
    arg='./ActuatorDef/Act_ActuatorDef.dat'
    fid=10
    OPEN(fid,file=trim(arg),status="old",action="read", iostat=io)

    READ(fid,'(A132)') char1
    READ(fid,'(A132)') char1
    READ(fid,'(A132)') char1
    READ(fid,'(A132)') char1
    READ(fid,*) farm1%num_turbines

if(t==0.d0 .and. OMP_GET_THREAD_NUM()==0)then
    print *,'-------------------------'
    print *,achar(27)//'[34m READING FARM DEFINITION'
    write(*,*) "Number of turbines in farm:", farm1%num_turbines
endif

    READ(fid,'(A132)') char1


    allocate(farm1%turbine(farm1%num_turbines))

     do i = 1, farm1%num_turbines
        READ(fid,*) farm1%turbine(i)%hub_cood_x, farm1%turbine(i)%hub_cood_y, farm1%turbine(i)%hub_cood_z
     ENDDO

     READ(fid,'(A132)') char1

     do i = 1, farm1%num_turbines
        READ(fid,*) farm1%turbine(i)%radius
     ENDDO

     READ(fid,'(A132)') char1

     do i = 1, farm1%num_turbines
        READ(fid,*) farm1%turbine(i)%normal_x, farm1%turbine(i)%normal_y, farm1%turbine(i)%normal_z
     ENDDO
    
    ! write(*,*) normal_x(:), normal_y(:), normal_z(:)
     READ(fid,'(A132)') char1

     do i = 1, farm1%num_turbines
        READ(fid,*) farm1%turbine(i)%rot_speed
     ENDDO

     READ(fid,'(A132)') char1

     do i = 1, farm1%num_turbines
        READ(fid,*) farm1%turbine(i)%blade_pitch
     ENDDO
    

     READ(fid,'(A132)') char1
     READ(fid,'(A132)') char1

     ! Read blade info, we assume all 3 blades are the same for one turbine
     READ(fid,*) farm1%turbine(1)%num_blade_sections

if(t==0.d0 .and. OMP_GET_THREAD_NUM()==0)then
    write(*,*) "Number of blade sections:", farm1%turbine(1)%num_blade_sections
endif
     associate (num_blade_sections => farm1%turbine(1)%num_blade_sections)

     READ(fid,'(A132)') char1

     do i=1, farm1%num_turbines
      do j=1, farm1%turbine(i)%num_blades
     allocate( farm1%turbine(i)%blade(j)%r_R(num_blade_sections),farm1%turbine(i)%blade(j)%chord(num_blade_sections), &
     farm1%turbine(i)%blade(j)%twist(num_blade_sections), farm1%turbine(i)%blade(j)%azimuth_angle(3), &
     farm1%turbine(i)%blade(j)%num_airfoils(num_blade_sections), &
     farm1%turbine(i)%blade(j)%airfoil_files(num_blade_sections,5),farm1%turbine(i)%blade(j)%airfoil(num_blade_sections), &
     farm1%turbine(i)%blade(j)%local_velocity(num_blade_sections), farm1%turbine(i)%blade(j)%local_angle(num_blade_sections), &
     farm1%turbine(i)%blade(j)%local_lift(num_blade_sections), farm1%turbine(i)%blade(j)%local_drag(num_blade_sections), &
     farm1%turbine(i)%blade(j)%point_xyz_loc(num_blade_sections,3),farm1%turbine(i)%blade(j)%local_torque(num_blade_sections), &
     farm1%turbine(i)%blade(j)%local_thrust(num_blade_sections),farm1%turbine(i)%blade(j)%local_root_bending(num_blade_sections)) 
     ! max 5 airfoils file names per section

         do k=1, num_blade_sections
            farm1%turbine(i)%blade(j)%airfoil_files(k,:)=' '
         enddo
      ENDDO
   enddo

    endassociate

   do i = 1, farm1%turbine(1)%num_blade_sections
      READ(fid,*) farm1%turbine(1)%blade(1)%r_R(i), farm1%turbine(1)%blade(1)%chord(i), &
                  farm1%turbine(1)%blade(1)%twist(i), farm1%turbine(1)%blade(1)%num_airfoils(i)
                    
      do j = 1, farm1%turbine(1)%blade(1)%num_airfoils(i)   
            READ(fid,*) farm1%turbine(1)%blade(1)%airfoil_files(i,j)  
      enddo

      ! azimuthal angle for the 3 blades
      farm1%turbine(1)%blade(1)%azimuth_angle(1)=0.d0
      farm1%turbine(1)%blade(1)%azimuth_angle(2)=pi*2.d0/3.
      farm1%turbine(1)%blade(1)%azimuth_angle(3)=pi*4.d0/3.
   ENDDO
     
     ! all turbines have the same blades
     !do i=1, farm1%num_turbines
     !    do j=1, farm1%turbine(i)%num_blades
     !       farm1%turbine(i)%blade(j)=farm1%turbine(1)%blade(1)
     !    ENDDO
     ! enddo
     ! write(*,*) "All turbines have the same blades"
   !  write(*,*) farm1%turbine(1)%blade(1)%airfoil_files(2,2)

! read numerical parameters
     READ(fid,'(A132)') char1
     READ(fid,'(A132)') char1
     READ(fid,'(A132)') char1
     READ(fid,'(A132)') char1

     READ(fid,*) epsil    

     READ(fid,'(A132)') char1
     READ(fid,*) c1,c2

if(t==0.d0 .and. OMP_GET_THREAD_NUM()==0)then
	print*,'Gaussian value for actuator line',epsil
	print*,'Tip correction constants', c1,c2
endif

    close(fid)
    !print *,'-------------------------'

    do i = 1, farm1%turbine(1)%num_blade_sections

      arg=trim('./ActuatorDef/'//trim(farm1%turbine(1)%blade(1)%airfoil_files(i,1)))
      !print*, 'reading: ', trim(arg)

      fid=10
      OPEN(fid,file=trim(arg),status="old",action="read", iostat=io)

      READ(fid,'(A132)') char1

      READ(fid,*) farm1%turbine(1)%blade(1)%airfoil(i)%num_aoa

if(t==0.d0 .and. OMP_GET_THREAD_NUM()==0)then
      print *,'-------------------------'
      print *,achar(27)//'[34m READING FARM AIRFOIL DATA (Cl-Cd)'
      print*, 'reading: ', trim(arg)
      write(*,*) 'The number of AoA in the file is: ', farm1%turbine(1)%blade(1)%airfoil(i)%num_aoa,' '//achar(27)//'[0m '	
endif
    
      READ(fid,'(A132)') char1

    associate (num_aoa => farm1%turbine(1)%blade(1)%airfoil(i)%num_aoa)

    do ii=1, farm1%num_turbines
      do j=1, farm1%turbine(ii)%num_blades
         !do k=1, farm1%turbine(1)%blade(1)%num_airfoils(j) ! this needs changing if many airfoils per radial section
            allocate( farm1%turbine(ii)%blade(j)%airfoil(i)%aoa(num_aoa), &
                  farm1%turbine(ii)%blade(j)%airfoil(i)%cl(num_aoa), &
                  farm1%turbine(ii)%blade(j)%airfoil(i)%cd(num_aoa))
         !enddo
      ENDDO
   enddo

   endassociate

    do ii = 1,  farm1%turbine(1)%blade(1)%airfoil(i)%num_aoa
         READ(fid,*) farm1%turbine(1)%blade(1)%airfoil(i)%aoa(ii), farm1%turbine(1)%blade(1)%airfoil(i)%cl(ii), &
                     farm1%turbine(1)%blade(1)%airfoil(i)%cd(ii)
    enddo

    !do ii = 1,  farm1%turbine(1)%blade(1)%airfoil(i)%num_aoa
    !     write(*,fmt='(F6.2,F6.2,F6.2)') farm1%turbine(1)%blade(1)%airfoil(i)%aoa(ii), farm1%turbine(1)%blade(1)%airfoil(i)%cl(ii), &
    !     farm1%turbine(1)%blade(1)%airfoil(i)%cd(ii)
    !enddo



    !all airfoils of all blades of all turbines are the same
    !do ii=1, farm1%num_turbines
    !  do j=1, farm1%turbine(i)%num_blades
    !     do k=1, farm1%turbine(1)%blade(1)%num_airfoils(j) 
    !     farm1%turbine(ii)%blade(j)%airfoil(k)=farm1%turbine(1)%blade(1)%airfoil(1)
    !     enddo
    !  ENDDO
    !enddo
      
    close(fid)
   enddo ! number of blade sections
!$omp end critical

        !all airfoils of all blades of all turbines are the same
   do ii=1, farm1%num_turbines
      do j=1, farm1%turbine(ii)%num_blades 
         farm1%turbine(ii)%blade(j)=farm1%turbine(1)%blade(1)
         enddo
   enddo


   Cl=0.0
   Cd=0.0 
   !wind_speed=1.0

   ! initial transcient (avoid sudden start)
   !t_init = 1. ! el tiempo cuando el funcionamiento es ya normal

   !if (t < t_init) then
   !    red_factor = 1 !t*t_init 
   !else
   !    red_factor = 1                
   !end if

   theta=farm1%turbine(1)%rot_speed*t

   !x coodrinate of every acutator line point
   farm1%turbine(1)%blade(3)%point_xyz_loc(:,1)=0.d0

   ! y,z coodrinate of every acutator line point

!$omp do schedule(runtime)private(i)
   do i = 1, farm1%turbine(1)%num_blade_sections

      farm1%turbine(1)%blade(1)%point_xyz_loc(i,2) =  farm1%turbine(1)%hub_cood_y+farm1%turbine(1)%radius*real(i) *cos(theta) &
                                                      /farm1%turbine(1)%num_blade_sections
      farm1%turbine(1)%blade(1)%point_xyz_loc(i,3) =   farm1%turbine(1)%hub_cood_z+farm1%turbine(1)%radius*real(i) *sin(theta) &
                                                      /farm1%turbine(1)%num_blade_sections

      farm1%turbine(1)%blade(2)%point_xyz_loc(i,2) =   farm1%turbine(1)%hub_cood_y+farm1%turbine(1)%radius*real(i) *cos(theta+2.0*PI/3.0) &
                                                      /farm1%turbine(1)%num_blade_sections
      farm1%turbine(1)%blade(2)%point_xyz_loc(i,3) =   farm1%turbine(1)%hub_cood_z+farm1%turbine(1)%radius*real(i) *sin(theta+2.0*PI/3.0) &
                                                      /farm1%turbine(1)%num_blade_sections

      farm1%turbine(1)%blade(3)%point_xyz_loc(i,2) =   farm1%turbine(1)%hub_cood_y+farm1%turbine(1)%radius*real(i) *cos(theta+4.0*PI/3.0) &
                                                      /farm1%turbine(1)%num_blade_sections
      farm1%turbine(1)%blade(3)%point_xyz_loc(i,3) =   farm1%turbine(1)%hub_cood_z+farm1%turbine(1)%radius*real(i) *sin(theta+4.0*PI/3.0) &
                                                      /farm1%turbine(1)%num_blade_sections
   enddo


!        *****************************************************************************************************************************
!        Compute contributions to source term
!        ATTENTION: This is deactivated for child multigrid meshes since they have specially tailored source terms (already computed).
!                   If you are going to add contributions to the source term, do it adding to e % storage % S_NS inside the condition!
!        *****************************************************************************************************************************
!
        if (.not. mesh % child) then
!
!
!           Add physical source term
!           ************************
!$omp do schedule(runtime) private(eID,i,j,k,ii,jj)
            do eID = 1, mesh % no_of_elements
               associate ( e => mesh % elements(eID) )
               do k = 0, e % Nxyz(3)   ; do j = 0, e % Nxyz(2) ; do i = 0, e % Nxyz(1)

! turbine is pointing backwards as x poistive
if((e % geom % x(2,i,j,k)-farm1%turbine(1)%hub_cood_y)**2+(e % geom % x(3,i,j,k)-farm1%turbine(1)%hub_cood_z)**2<=farm1%turbine(1)%radius .and. (e % geom % x(1,i,j,k)<0.2*farm1%turbine(1)%radius .and. e %geom % x(1,i,j,k)>-0.2*farm1%turbine(1)%radius)) then
! 20% of the radius for max of the rotor thikness

farm1%turbine(1)%blade_torque(:)=0.d0
farm1%turbine(1)%blade_thrust(:)=0.d0
farm1%turbine(1)%blade_root_bending(:)=0.d0


   do jj = 1, farm1%turbine(1)%num_blades
      
      do ii = 1, farm1%turbine(1)%num_blade_sections

   	wind_speed=e % storage %Q(2,i,j,k)/e % storage %Q(1,i,j,k) ! wind goes in the x-direction
	!wind_speed=1.d0
	
	density=e % storage %Q(1,i,j,k)

         farm1%turbine(1)%blade(jj)%local_velocity(ii)= &
               sqrt(farm1%turbine(1)%rot_speed*farm1%turbine(1)%blade(jj)%r_R(ii)**2 &
               +wind_speed**2)

         farm1%turbine(1)%blade(jj)%local_angle(ii)=atan(wind_speed/farm1%turbine(1)%blade(jj)%local_velocity(ii)) & 
                                                -farm1%turbine(1)%blade(jj)%twist(ii) &
                                                -farm1%turbine(1)%blade_pitch

!tip correction
g1_func=exp(-c1*(farm1%turbine(1)%num_blades*farm1%turbine(1)%blade(jj)%local_velocity(ii)/wind_speed-c2))+0.1
tip_correct=2.d0/pi*acos(exp(-g1_func*farm1%turbine(1)%num_blades*(farm1%turbine(1)%radius-farm1%turbine(1)%blade(jj)%r_R(ii))/(2d0*farm1%turbine(1)%blade(jj)%r_R(ii)*sin(farm1%turbine(1)%blade(jj)%local_angle(ii)))))


         ! ITERPOLATE     
         call Get_Cl_Cl_from_airfoil_data(farm1%turbine(1)%blade(jj)%airfoil(ii), aoa, Cl, Cd)

         
         interp = exp(-(norm2([e % geom %x(1,i,j,k),e % geom %x(2,i,j,k),e % geom %x(3,i,j,k)]-farm1%turbine(1)%blade(jj)%point_xyz_loc(ii,:))**2) & 
                  /epsil/epsil)/sqrt((2.*pi)**3)/epsil**3

         farm1%turbine(1)%blade(jj)%local_lift(ii)=-farm1%turbine(1)%blade(jj)%local_velocity(ii)**2 &
                                                         *0.5*density*Cl*interp

         farm1%turbine(1)%blade(jj)%local_drag(ii)=-farm1%turbine(1)%blade(jj)%local_velocity(ii)**2 &
                                                         *0.5*density*Cd*interp

         Rotor_force = + farm1%turbine(1)%blade(jj)%local_lift(ii)*sin(farm1%turbine(1)%blade(jj)%local_angle(ii)) &
                        - farm1%turbine(1)%blade(jj)%local_drag(ii)*cos(farm1%turbine(1)%blade(jj)%local_angle(ii))

                                    
! the rotor normal is poitning in the x direction
              
         e % storage % S_NS(2,i,j,k)= tip_correct*farm1%turbine(1)%blade(jj)%local_lift(ii)*cos(farm1%turbine(1)%blade(jj)%local_angle(ii)) + farm1%turbine(1)%blade(jj)%local_drag(ii)*sin(farm1%turbine(1)%blade(jj)%local_angle(ii))

	e % storage % S_NS(3,i,j,k)=- tip_correct*Rotor_force*cos(farm1%turbine(1)%rot_speed*t+farm1%turbine(1)%blade(jj)%azimuth_angle(jj))

	e % storage % S_NS(4,i,j,k)=- tip_correct*Rotor_force*sin(farm1%turbine(1)%rot_speed*t+farm1%turbine(1)%blade(jj)%azimuth_angle(jj))


farm1%turbine(1)%blade(jj)%local_thrust(ii)=e % storage % S_NS(2,i,j,k)
farm1%turbine(1)%blade(jj)%local_torque(ii)=sqrt(e % storage % S_NS(3,i,j,k)**2+e % storage % S_NS(4,i,j,k)**2)*farm1%turbine(1)%blade(jj)%r_R(ii)

farm1%turbine(1)%blade(jj)%local_root_bending(ii)=farm1%turbine(1)%blade(jj)%local_thrust(ii)*farm1%turbine(1)%blade(jj)%r_R(ii)

!print*, farm1%turbine(1)%blade(jj)%local_thrust(ii)
       enddo

farm1%turbine(1)%blade_thrust(jj)=sum(farm1%turbine(1)%blade(jj)%local_thrust)
farm1%turbine(1)%blade_torque(jj)=sum(farm1%turbine(1)%blade(jj)%local_torque)
farm1%turbine(1)%blade_root_bending(jj)=sum(farm1%turbine(1)%blade(jj)%local_root_bending)
   enddo

endif  ! coordenates

                  call UserDefinedSourceTermNS(e % geom % x(:,i,j,k), e % storage % Q(:,i,j,k), t, e % storage % S_NS(:,i,j,k), thermodynamics, dimensionless, refValues)


               end do                  ; end do                ; end do
               end associate
            end do
!$omp end do

!write output tirque thrust to file

if(OMP_GET_THREAD_NUM()==0)then
      fid=11
      arg='./RESULTS/Output_Actuator_Forces.dat'
  if(t==0.d0)then
      OPEN(fid,file=trim(arg), status="new",iostat=io)
      write(fid,*) 'time, thrust_1, blade_torque_1, blade_root_bending_1,thrust_2, blade_torque_12 blade_root_bending_2,thrust_3, blade_torque_3, blade_root_bending_3'
      write(fid,"(10(2X,E11.4))") t, farm1%turbine(1)%blade_thrust(1),farm1%turbine(1)%blade_torque(1),farm1%turbine(1)%blade_root_bending(1), &
	farm1%turbine(1)%blade_thrust(2),farm1%turbine(1)%blade_torque(2),farm1%turbine(1)%blade_root_bending(2),&
	farm1%turbine(1)%blade_thrust(3),farm1%turbine(1)%blade_torque(3),farm1%turbine(1)%blade_root_bending(3)
      close(fid)
  else
      OPEN(fid,file=trim(arg), status="old",POSITION='APPEND',iostat=io)
      write(fid,"(10(2X,E11.4))") t, farm1%turbine(1)%blade_thrust(1),farm1%turbine(1)%blade_torque(1),farm1%turbine(1)%blade_root_bending(1), &
	farm1%turbine(1)%blade_thrust(2),farm1%turbine(1)%blade_torque(2),farm1%turbine(1)%blade_root_bending(2),&
	farm1%turbine(1)%blade_thrust(3),farm1%turbine(1)%blade_torque(3),farm1%turbine(1)%blade_root_bending(3)
      close(fid)
  endif
endif


!           Add Particles source
!           ********************
            if ( particles % active ) then            
!$omp do schedule(runtime)            
               do eID = 1, mesh % no_of_elements
                  associate ( e => mesh % elements(eID) )            
                     e % storage % S_NSP = 0.0_RP
                  end associate
               enddo 
!$omp end do             
               
!$omp do schedule(runtime) 
               do i = 1, particles % injection % injected + 1
                  if (particles % particle(i) % active) then 

                     associate ( eID => particles % particle(i) % eID )
                        
                     call particles % AddSource(i, mesh % elements( eID ), &
                        t, thermodynamics, dimensionless, refValues)

                     ! If this is uncommented, private(j) should be added to openmp.
                        !this commented section permits the computation of source term in neighbor elements
                     !do j = 1, 6
                     !   if (particles % particle(i) % mesh%elements( eID )%NumberOfConnections(j) > 0) then  
                     !      call particles % AddSource(i, &
                     !      mesh % elements( particles % particle(i) % mesh%elements( eID )%Connection(j)%ElementIDs(1)  ), &
                     !      t, thermodynamics, dimensionless, refValues)
                     !   else 
                     !      !
                     !   end if
                     !end do   
                     end associate   
                  endif 
               end do
!$omp end do
!$omp do schedule(runtime)  
               do eID = 1, mesh % no_of_elements
                  associate ( e => mesh % elements(eID) )            
                     e % storage % S_NS = e % storage % S_NS + e % storage % S_NSP
                  end associate
               enddo
!$omp end do
            end if
         end if !(.not. mesh % child)
!
!        ***********************
!        Now add the source term
!        ***********************
!$omp do schedule(runtime) private(i,j,k)
         do eID = 1, mesh % no_of_elements
            associate ( e => mesh % elements(eID) )
            do k = 0, e % Nxyz(3)   ; do j = 0, e % Nxyz(2) ; do i = 0, e % Nxyz(1)
               e % storage % QDot(:,i,j,k) = e % storage % QDot(:,i,j,k) + e % storage % S_NS(:,i,j,k)
            end do                  ; end do                ; end do
            end associate
         end do
!$omp end do


deallocate(farm1%turbine)

      end subroutine TimeDerivative_ComputeQDot
   
      subroutine compute_viscosity_at_faces(no_of_faces, no_of_sides, face_ids, mesh)
         implicit none
         integer, intent(in)           :: no_of_faces
         integer, intent(in)           :: no_of_sides
         integer, intent(in)           :: face_ids(no_of_faces)
         class(HexMesh), intent(inout) :: mesh
!
!        ---------------
!        Local variables
!        ---------------
!
         integer       :: iFace, i, j, side
         real(kind=RP) :: delta, mu_smag

         if (flowIsNavierStokes) then
!$omp do schedule(runtime) private(i,j)
            do iFace = 1, no_of_faces
               associate(f => mesh % faces(face_ids(iFace)))
               do j = 0, f % Nf(2) ; do i = 0, f % Nf(1)
                  do side = 1, no_of_sides
                      call get_laminar_mu_kappa(f % storage(side) % Q(:,i,j), f % storage(side) % mu_NS(1,i,j), f % storage(side) % mu_NS(2,i,j))
                  end do
               end do              ; end do         
               end associate
            end do
!$omp end do
         end if

         if ( LESModel % Active ) then
!$omp do schedule(runtime) private(i,j,delta,mu_smag)
            do iFace = 1, no_of_faces
               associate(f => mesh % faces(face_ids(iFace)))

               delta = sqrt(f % geom % surface / product(f % Nf + 1))
               do j = 0, f % Nf(2) ; do i = 0, f % Nf(1)
                  do side = 1, no_of_sides
                     call LESModel % ComputeViscosity(delta, f % geom % dWall(i,j), f % storage(side) % Q(:,i,j),   &
                                                                                    f % storage(side) % U_x(:,i,j), &
                                                                                    f % storage(side) % U_y(:,i,j), &
                                                                                    f % storage(side) % U_z(:,i,j), &
                                                                                    mu_smag)
                     f % storage(side) % mu_NS(1,i,j) = f % storage(side) % mu_NS(1,i,j) + mu_smag
                     f % storage(side) % mu_NS(2,i,j) = f % storage(side) % mu_NS(2,i,j) + mu_smag * dimensionless % mu_to_kappa
                  end do
               end do              ; end do         
               end associate
            end do
!$omp end do
         end if



      end subroutine compute_viscosity_at_faces
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
         procedure(UserDefinedSourceTermNS_f) :: UserDefinedSourceTermNS
!
!        ****************
!        Volume integrals
!        ****************
!
!$omp do schedule(runtime) 
         do eID = 1 , size(mesh % elements)
            call TimeDerivative_StrongVolumetricContribution(mesh, mesh % elements(eID) , t)
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
         
!$omp do schedule(runtime) private(i,j,k)
         do eID = 1, mesh % no_of_elements
            associate ( e => mesh % elements(eID) )
            do k = 0, e % Nxyz(3)   ; do j = 0, e % Nxyz(2) ; do i = 0, e % Nxyz(1)
               e % storage % QDot(:,i,j,k) = e % storage % QDot(:,i,j,k) + e % storage % S_NS(:,i,j,k)
            end do                  ; end do                ; end do
            end associate
         end do
!$omp end do
         
      end subroutine TimeDerivative_ComputeQDotIsolated
!
!///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
!
      subroutine TimeDerivative_StrongVolumetricContribution(mesh, e, t)
         use HexMeshClass
         use ElementClass
         implicit none
         type(HexMesh)             :: mesh
         type(Element)             :: e
         real(kind=RP), intent(in) :: t

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
         call HyperbolicDiscretization % ComputeInnerFluxes ( e , EulerFlux, inviscidContravariantFlux ) 
!
!        Compute viscous contravariant flux
!        ----------------------------------
         call ViscousDiscretization  % ComputeInnerFluxes ( NCONS, NGRAD, ViscousFlux, GetNSViscosity, e, viscousContravariantFlux) 
!
!        Compute the SVV dissipation
!        ---------------------------
         if ( .not. SVV % enabled ) then
            SVVcontravariantFlux = 0.0_RP
         else
            call SVV % ComputeInnerFluxes(mesh, e, SVVContravariantFlux)
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
            e % storage % QDot = ScalarStrongIntegrals % StdVolumeGreen ( e , NCONS, contravariantFlux ) 

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
      subroutine TimeDerivative_VolumetricContribution(mesh, e, t)
         use HexMeshClass
         use ElementClass
         implicit none
         type(HexMesh)             :: mesh
         type(Element)             :: e
         real(kind=RP), intent(in) :: t

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
         call HyperbolicDiscretization % ComputeInnerFluxes ( e , EulerFlux, inviscidContravariantFlux ) 
!
!        Compute viscous contravariant flux
!        ----------------------------------
         if (flowIsNavierStokes) then
            call ViscousDiscretization  % ComputeInnerFluxes ( NCONS, NGRAD, ViscousFlux, GetNSViscosity, e , viscousContravariantFlux) 
!   
!           Compute the SVV dissipation
!           ---------------------------
            if ( .not. SVV % enabled ) then
               SVVcontravariantFlux = 0.0_RP
            else
               call SVV % ComputeInnerFluxes(mesh, e, SVVContravariantFlux)
            end if
         else
            viscousContravariantFlux = 0.0_RP
            SVVcontravariantFlux = 0.0_RP
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
         real(kind=RP) :: mu_left(3), mu_right(3)
        integer        :: Sidearray(2)
         
         visc_flux = 0._RP
         
         if (flowIsNavierStokes) then
            DO j = 0, f % Nf(2)
               DO i = 0, f % Nf(1)

                  mu_left(1) = f % storage(1) % mu_NS(1,i,j)
                  mu_left(2) = 0.0_RP
                  mu_left(3) = f % storage(1) % mu_NS(2,i,j)

                  mu_right(1) = f % storage(2) % mu_NS(1,i,j)
                  mu_right(2) = 0.0_RP
                  mu_right(3) = f % storage(2) % mu_NS(2,i,j)
!      
!                 --------------
!                 Viscous fluxes
!                 --------------
!      
                  CALL ViscousDiscretization % RiemannSolver(nEqn = NCONS, nGradEqn = NGRAD, &
                                                     EllipticFlux = ViscousFlux, &
                                                     f = f, &
                                                     QLeft = f % storage(1) % Q(:,i,j), &
                                                     QRight = f % storage(2) % Q(:,i,j), &
                                                     U_xLeft = f % storage(1) % U_x(:,i,j), &
                                                     U_yLeft = f % storage(1) % U_y(:,i,j), &
                                                     U_zLeft = f % storage(1) % U_z(:,i,j), &
                                                     U_xRight = f % storage(2) % U_x(:,i,j), &
                                                     U_yRight = f % storage(2) % U_y(:,i,j), &
                                                     U_zRight = f % storage(2) % U_z(:,i,j), &
                                                     mu_left = mu_left, mu_right = mu_right, &
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
         Sidearray = (/1,2/)
         call f % ProjectFluxToElements(NCONS, flux, Sidearray)

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
         real(kind=RP) :: mu_left(3), mu_right(3)
         integer       :: Sidearray(2)
         
         visc_flux = 0._RP
         if (flowIsNavierStokes) then
            DO j = 0, f % Nf(2)
               DO i = 0, f % Nf(1)

                  mu_left(1) = f % storage(1) % mu_NS(1,i,j)
                  mu_left(2) = 0.0_RP
                  mu_left(3) = f % storage(1) % mu_NS(2,i,j)

                  mu_right(1) = f % storage(2) % mu_NS(1,i,j)
                  mu_right(2) = 0.0_RP
                  mu_right(3) = f % storage(2) % mu_NS(2,i,j)
!      
!                 --------------
!                 Viscous fluxes
!                 --------------
!      
                  CALL ViscousDiscretization % RiemannSolver(nEqn = NCONS, nGradEqn = NGRAD, &
                                                     EllipticFlux = ViscousFlux, &
                                                     f = f, &
                                                     QLeft = f % storage(1) % Q(:,i,j), &
                                                     QRight = f % storage(2) % Q(:,i,j), &
                                                     U_xLeft = f % storage(1) % U_x(:,i,j), &
                                                     U_yLeft = f % storage(1) % U_y(:,i,j), &
                                                     U_zLeft = f % storage(1) % U_z(:,i,j), &
                                                     U_xRight = f % storage(2) % U_x(:,i,j), &
                                                     U_yRight = f % storage(2) % U_y(:,i,j), &
                                                     U_zRight = f % storage(2) % U_z(:,i,j), &
                                                     mu_left  = mu_left, &
                                                     mu_right = mu_right, &
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
!        Return the flux to elements: The sign in eR % storage % FstarB has already been accouted.
!        ---------------------------
!
         thisSide = maxloc(f % elementIDs, dim = 1)
         
         Sidearray = (/thisSide, HMESH_NONE/)
         call f % ProjectFluxToElements(NCONS, flux, Sidearray )

      end subroutine ComputeMPIFaceFlux_NS

      SUBROUTINE computeBoundaryFlux_NS(f, time, mesh)
      USE ElementClass
      use FaceClass
      USE RiemannSolvers_NS
      use WallFunctionBC
      use WallFunctionConnectivity
      IMPLICIT NONE
!
!     ---------
!     Arguments
!     ---------
!
      type(Face),    intent(inout) :: f
      REAL(KIND=RP)                :: time
      type(HexMesh), intent(in)    :: mesh
!
!     ---------------
!     Local variables
!     ---------------
!
      INTEGER                         :: i, j
      INTEGER, DIMENSION(2)           :: N
      REAL(KIND=RP)                   :: inv_flux(NCONS)
      real(kind=RP)                   :: visc_flux(NCONS, 0:f % Nf(1), 0:f % Nf(2))
      real(kind=RP)                   :: fStar(NCONS, 0:f % Nf(1), 0: f % Nf(2))
      real(kind=RP)                   :: mu, kappa, beta, delta
      real(kind=RP)                   :: fv_3d(NCONS,NDIM)
      integer                         :: Sidearray(2)
      logical                         :: useWallFuncFace
      real(kind=RP)                   :: wallFunV(NDIM, 0:f % Nf(1), 0:f % Nf(2))
      real(kind=RP)                   :: wallFunRho(0:f % Nf(1), 0:f % Nf(2))
      real(kind=RP)                   :: wallFunMu(0:f % Nf(1), 0:f % Nf(2))
      real(kind=RP)                   :: wallFunY(0:f % Nf(1), 0:f % Nf(2))
!
!     -------------------
!     Get external states
!     -------------------
!
      do j = 0, f % Nf(2)  ; do i = 0, f % Nf(1)
         f % storage(2) % Q(:,i,j) = f % storage(1) % Q(:,i,j)
         CALL BCs(f % zone) % bc % FlowState(f % geom % x(:,i,j), &
                                      time, &
                                      f % geom % normal(:,i,j), &
                                      f % storage(2) % Q(:,i,j))

      end do               ; end do

      if ( flowIsNavierStokes ) then

          useWallFuncFace = .false.
          if (useWallFunc) then
              do i = 1, len(wallFunBCs)
                  if (trim(wallFunBCs(i)) .eq. trim(f % boundaryName)) then
                      useWallFuncFace = .true.
                      exit
                  end if
              end do
          end if
          if (useWallFuncFace) then
              call WallFunctionGatherFlowVariables(mesh, f, wallFunV, wallFunRho, wallFunMu, wallFunY)
          end if

         DO j = 0, f % Nf(2)
            DO i = 0, f % Nf(1)
               mu    = f % storage(1) % mu_ns(1,i,j)
               beta  = 0.0_RP
               kappa = f % storage(1) % mu_ns(2,i,j)

               call ViscousFlux(NCONS,NGRAD,f % storage(1) % Q(:,i,j), &
                                            f % storage(1) % U_x(:,i,j), &
                                            f % storage(1) % U_y(:,i,j), &
                                            f % storage(1) % U_z(:,i,j), &
                                            mu, beta, kappa, fv_3d)

               visc_flux(:,i,j) =   fv_3d(:,IX)*f % geom % normal(IX,i,j) &
                                  + fv_3d(:,IY)*f % geom % normal(IY,i,j) &
                                  + fv_3d(:,IZ)*f % geom % normal(IZ,i,j) 
            
               if (useWallFuncFace) then
                   call WallViscousFlux(wallFunV(:,i,j), wallFunY(i,j), f % geom % normal(:,i,j), wallFunRho(i,j), wallFunMu(i,j), visc_flux(:,i,j))
               end if 

               CALL BCs(f % zone) % bc % FlowNeumann(&
                                              f % geom % x(:,i,j), &
                                              time, &
                                              f % geom % normal(:,i,j), &
                                              f % storage(1) % Q(:,i,j), &
                                              f % storage(1) % U_x(:,i,j), &
                                              f % storage(1) % U_y(:,i,j), &
                                              f % storage(1) % U_z(:,i,j), &
                                              visc_flux(:,i,j))

            end do
         end do
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

      Sidearray = (/1, HMESH_NONE/)
      call f % ProjectFluxToElements(NCONS, fStar, Sidearray)

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
         real(kind=RP) :: mu_left(3), mu_right(3)
         integer       :: Sidearray(2)
!
!        ----------
!        SVV fluxes
!        ----------
!
         if ( SVV % enabled ) then 
            SVV_Flux = 0.5_RP * (f % storage(1) % HFlux + f % storage(2) % HFlux)
         else
            SVV_Flux = 0.0_RP
         end if

         DO j = 0, f % Nf(2)
            DO i = 0, f % Nf(1)

!      
!              --------------
!              Viscous fluxes
!              --------------
!      
               mu_left(1) = f % storage(1) % mu_ns(1,i,j)
               mu_left(2) = 0.0_RP
               mu_left(3) = f % storage(1) % mu_ns(2,i,j)

               mu_right(1) = f % storage(2) % mu_ns(1,i,j)
               mu_right(2) = 0.0_RP
               mu_right(3) = f % storage(2) % mu_ns(2,i,j)

               CALL ViscousDiscretization % RiemannSolver(nEqn = NCONS, nGradEqn = NGRAD, &
                                                  EllipticFlux = ViscousFlux, &
                                                  f = f, &
                                                  QLeft = f % storage(1) % Q(:,i,j), &
                                                  QRight = f % storage(2) % Q(:,i,j), &
                                                  U_xLeft = f % storage(1) % U_x(:,i,j), &
                                                  U_yLeft = f % storage(1) % U_y(:,i,j), &
                                                  U_zLeft = f % storage(1) % U_z(:,i,j), &
                                                  U_xRight = f % storage(2) % U_x(:,i,j), &
                                                  U_yRight = f % storage(2) % U_y(:,i,j), &
                                                  U_zRight = f % storage(2) % U_z(:,i,j), &
                                                  mu_left  = mu_left, &
                                                  mu_right = mu_right, &
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
               flux(:,i,j) = (inv_flux(:,i,j) - visc_flux(:,i,j)) * f % geom % jacobian(i,j) - SVV_flux(:,i,j)
               
            END DO   
         END DO  
!
!        ---------------------------
!        Return the flux to elements
!        ---------------------------
!
         Sidearray = (/1,2/)
         call f % ProjectFluxToElements(NCONS, flux, Sidearray )

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
         integer       :: Sidearray(2)
         real(kind=RP) :: mu(3)
!
!        ----------
!        SVV fluxes
!        ----------
!
         if ( SVV % enabled ) then
            print*, "SVV Not configured with MPI"
            errorMessage(STD_OUT)
            stop
         end if

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
               mu(1) = f % storage(1) % mu_ns(1,i,j)
               mu(2) = 0.0_RP 
               mu(3) = f % storage(1) % mu_ns(2,i,j)

               CALL ViscousDiscretization % RiemannSolver(nEqn = NCONS, nGradEqn = NGRAD, &
                                                  EllipticFlux = ViscousFlux, &
                                                  f = f, &
                                                  QLeft = f % storage(1) % Q(:,i,j), &
                                                  QRight = f % storage(2) % Q(:,i,j), &
                                                  U_xLeft = f % storage(1) % U_x(:,i,j), &
                                                  U_yLeft = f % storage(1) % U_y(:,i,j), &
                                                  U_zLeft = f % storage(1) % U_z(:,i,j), &
                                                  U_xRight = f % storage(2) % U_x(:,i,j), &
                                                  U_yRight = f % storage(2) % U_y(:,i,j), &
                                                  U_zRight = f % storage(2) % U_z(:,i,j), &
                                                  mu_left  = mu, &
                                                  mu_right = mu, &
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
         
         Sidearray = (/thisSide, HMESH_NONE/)
         call f % ProjectFluxToElements(NCONS, flux, Sidearray)

      end subroutine ComputeMPIFaceFlux_SVV

      SUBROUTINE computeBoundaryFlux_SVV(f, time, mesh)
      USE ElementClass
      use FaceClass
      USE RiemannSolvers_NS
      IMPLICIT NONE
!
!     ---------
!     Arguments
!     ---------
!
      type(Face),    intent(inout) :: f
      REAL(KIND=RP)                :: time
        type(HexMesh), intent(in)  :: mesh
!
!     ---------------
!     Local variables
!     ---------------
!
      INTEGER               :: i, j
      INTEGER, DIMENSION(2) :: N
      REAL(KIND=RP)         :: inv_flux(NCONS)
      real(kind=RP)         :: visc_flux(NCONS, 0:f % Nf(1), 0:f % Nf(2))
      real(kind=RP)         :: SVV_flux(NCONS, 0:f % Nf(1), 0:f % Nf(2))
      real(kind=RP)         :: fStar(NCONS, 0:f % Nf(1), 0: f % Nf(2))
      real(kind=RP)         :: mu, beta, kappa, delta
      real(kind=RP)         :: tauSGS(NDIM,NDIM), qSGS(NDIM), fv_3d(NCONS,NDIM)
      integer               :: Sidearray(2)

      if ( SVV % enabled ) then
         do j = 0, f % Nf(2) ; do i = 0, f % Nf(1)
            SVV_flux(:,i,j) = f % storage(1) % Hflux(:,i,j) / f % geom % jacobian(i,j)
         end do              ; end do
      else
         SVV_flux = 0.0_RP
      end if

!
!     -------------------
!     Get external states
!     -------------------
!
      do j = 0, f % Nf(2)  ; do i = 0, f % Nf(1)
         f % storage(2) % Q(:,i,j) = f % storage(1) % Q(:,i,j)
         CALL BCs(f % zone) % bc % FlowState(f % geom % x(:,i,j), &
                                      time, &
                                      f % geom % normal(:,i,j), &
                                      f % storage(2) % Q(:,i,j))

      end do               ; end do

      if ( flowIsNavierStokes ) then
         DO j = 0, f % Nf(2)
            DO i = 0, f % Nf(1)
               mu    = f % storage(1) % mu_ns(1,i,j)
               beta  = 0.0_RP
               kappa = f % storage(1) % mu_ns(2,i,j)

               call ViscousFlux(NCONS,NGRAD,f % storage(1) % Q(:,i,j), &
                                            f % storage(1) % U_x(:,i,j), &
                                            f % storage(1) % U_y(:,i,j), &
                                            f % storage(1) % U_z(:,i,j), &
                                            mu, beta, kappa, fv_3d)

               visc_flux(:,i,j) =   fv_3d(:,IX)*f % geom % normal(IX,i,j) &
                                  + fv_3d(:,IY)*f % geom % normal(IY,i,j) &
                                  + fv_3d(:,IZ)*f % geom % normal(IZ,i,j) 

               visc_flux(:,i,j) = visc_flux(:,i,j) + SVV_flux(:,i,j)

               CALL BCs(f % zone) % bc % FlowNeumann(&
                                              f % geom % x(:,i,j), &
                                              time, &
                                              f % geom % normal(:,i,j), &
                                              f % storage(1) % Q(:,i,j), &
                                              f % storage(1) % U_x(:,i,j), &
                                              f % storage(1) % U_y(:,i,j), &
                                              f % storage(1) % U_z(:,i,j), &
                                              visc_flux(:,i,j))

            end do
         end do
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

      Sidearray = (/1, HMESH_NONE/)
      call f % ProjectFluxToElements(NCONS, fStar, Sidearray)

      END SUBROUTINE computeBoundaryFlux_SVV


! linear interpolation given two points; returns y for new_x following line coefs (a,b) with y=ax+b
function InterpolateAirfoilData(x1,x2,y1,y2,new_x)
   implicit none
    
   real(KIND=RP), intent(in)    :: x1, x2, y1, y2, new_x
   real(KIND=RP)                :: a, b, InterpolateAirfoilData

    if(abs(x1-x2)<1.0e-6) then 
      a=100.d0
   else
      a=(y1- y2)/(x1- x2)
   endif
    b= y1-a*x1;
    InterpolateAirfoilData=a*new_x+b
end function

subroutine Get_Cl_Cl_from_airfoil_data(airfoil1, aoa, Cl_out, Cd_out)
         implicit none
      

!         type airfoil
!         integer                  :: num_aoa 
!         real(KIND=RP), allocatable      :: aoa(:)  ! in rad
!         real(KIND=RP), allocatable      :: cl(:)  ! in rad
!         real(KIND=RP), allocatable      :: cd(:)  ! in rad
!         end type

         type (airfoil), intent(in)    :: airfoil1
         real(KIND=RP), intent(in)            :: aoa
         real(KIND=RP), intent(inout)         :: Cl_out, Cd_out
         integer                       :: i

!         INTERFACE 
!            FUNCTION InterpolateAirfoilData(x1,x2,y1,y2,new_x)
!               real(KIND=RP)                :: InterpolateAirfoilData
!               real(KIND=RP), intent(in)    :: x1, x2, y1, y2, new_x
!            END FUNCTION 
!         END INTERFACE

               
  do i=1, airfoil1%num_aoa-1
      if (airfoil1%aoa(i+1)>=aoa .and. airfoil1%aoa(i)<=aoa ) then
         Cl_out=InterpolateAirfoilData(airfoil1%aoa(i),airfoil1%aoa(i+1),airfoil1%cl(i),airfoil1%cl(i+1),aoa)
         Cd_out=InterpolateAirfoilData(airfoil1%aoa(i),airfoil1%aoa(i+1),airfoil1%cd(i),airfoil1%cd(i+1),aoa)
         exit
      else
         Cl_out=0.d0
         Cd_out=0.d0
      endif
  end do
end subroutine


end module SpatialDiscretization
