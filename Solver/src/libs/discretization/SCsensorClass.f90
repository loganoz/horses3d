#include "Includes.h"
module SCsensorClass

   use SMConstants,       only: RP, PI, STD_OUT, NDIM, IX, IY, IZ, LINE_LENGTH
   use DGSEMClass,        only: ComputeTimeDerivative_f, DGSem
   use PhysicsStorage,    only: grad_vars, GRADVARS_STATE, GRADVARS_ENTROPY, NCONS
   use NodalStorageClass, only: NodalStorage
   use Utilities,         only: toLower
   use Clustering,        only: GMM

   use ShockCapturingKeywords

   implicit none

   private

   public :: SCsensor_t
   public :: Set_SCsensor
   public :: Destruct_SCsensor

   type :: TruncationError_t
      integer                  :: Nmin       !< Min order N
      integer                  :: deltaN     !< N' of coarser mesh is max(Nmin, N-deltaN)
      integer                  :: derivType  !< Whether the time derivative isolated or not
      type(DGSem), allocatable :: coarseSem  !< DGSEM for the coarser mesh
      procedure(ComputeTimeDerivative_f), nopass, pointer :: TimeDerivative => null()
   end type TruncationError_t

   type :: SCsensor_t

         integer  :: sens_type  ! ID of the sensor type

         real(RP) :: s0         !< Centre of the sensor scaling
         real(RP) :: ds         !< Bandwith of the sensor scaling
         real(RP) :: ds2        !< Half-bandwidth
         real(RP) :: low        !< Lower threshold
         real(RP) :: high       !< Upper threshold
         integer  :: sVar       !< Variable used as input for the sensor
         integer  :: nclusters  !< Number of clusters (not for all the sensors)
         integer  :: min_steps  !< Minimum number of steps before the sensor is deactivated

         type(TruncationError_t), allocatable :: TEestim  !< Truncation error estimation

         procedure(Compute_Int), pointer :: Compute_Raw => null()

      contains

         procedure :: Describe => Describe_SCsensor
         procedure :: Compute  => Compute_SCsensor
         final     :: Destruct_SCsensor

   end type SCsensor_t
!
!  Interfaces
!  ----------
   abstract interface
      subroutine Compute_Int(this, sem, t)
         import RP, SCsensor_t, DGsem
         class(SCsensor_t), target, intent(inout) :: this
         type(DGSem),       target, intent(inout) :: sem
         real(RP),                  intent(in)    :: t
      end subroutine Compute_Int
   end interface
!
!  Private variables
!  -----------------
   integer, parameter :: ISOLATED_TE     = 0
   integer, parameter :: NON_ISOLATED_TE = 1
!
!  ========
   contains
!  ========
!
   subroutine Set_SCsensor(sensor, controlVariables, sem, minSteps, &
                           ComputeTimeDerivative, ComputeTimeDerivativeIsolated)
!
!     -------
!     Modules
!     -------
      use FTValueDictionaryClass
!
!     ---------
!     Interface
!     ---------
      implicit none
      type(SCsensor_t),                   intent(inout) :: sensor
      type(FTValueDictionary),            intent(in)    :: controlVariables
      type(DGSem),                        intent(in)    :: sem
      integer,                            intent(in)    :: minSteps
      procedure(ComputeTimeDerivative_f)                :: ComputeTimeDerivative
      procedure(ComputeTimeDerivative_f)                :: ComputeTimeDerivativeIsolated
!
!     ---------------
!     Local variables
!     ---------------
      character(len=:), allocatable :: sensorType
      character(len=:), allocatable :: sVar


      sensor % min_steps = minSteps
!
!     Sensor type
!     -----------
      if (controlVariables % containsKey(SC_SENSOR_KEY)) then
         sensorType = controlVariables % stringValueForKey(SC_SENSOR_KEY, LINE_LENGTH)
      else
         sensorType = SC_INTEGRAL_VAL
      end if
      call toLower(sensorType)

      select case (sensorType)
      case (SC_ZERO_VAL)
         sensor % sens_type = SC_ZERO_ID
         sensor % Compute_Raw => Sensor_zero

      case (SC_ONE_VAL)
         sensor % sens_type = SC_ONE_ID
         sensor % Compute_Raw => Sensor_one

      case (SC_INTEGRAL_VAL)
         sensor % sens_type = SC_INTEGRAL_ID
         sensor % Compute_Raw => Sensor_integral

      case (SC_MODAL_VAL)
         sensor % sens_type = SC_MODAL_ID
         sensor % Compute_Raw => Sensor_modal

      case (SC_TE_VAL)
         sensor % sens_type = SC_TE_ID
         sensor % Compute_Raw => Sensor_truncation
         call Construct_TEsensor(sensor, controlVariables, sem, &
                                 ComputeTimeDerivative, ComputeTimeDerivativeIsolated)

      case (SC_ALIAS_VAL)
         sensor % sens_type = SC_ALIAS_ID
         sensor % Compute_Raw => Sensor_aliasing

      case (SC_GMM_VAL)
         sensor % sens_type = SC_GMM_ID
         sensor % Compute_Raw => Sensor_GMM

      case default
         write(STD_OUT,*) "ERROR. The sensor type is unknown. Options are:"
         write(STD_OUT,*) '   * ', SC_ZERO_VAL
         write(STD_OUT,*) '   * ', SC_ONE_VAL
         write(STD_OUT,*) '   * ', SC_INTEGRAL_VAL
         write(STD_OUT,*) '   * ', SC_MODAL_VAL
         write(STD_OUT,*) '   * ', SC_TE_VAL
         write(STD_OUT,*) '   * ', SC_ALIAS_VAL
         write(STD_OUT,*) '   * ', SC_GMM_VAL
         stop

      end select
!
!     Options for the clustering sensor
!     ---------------------------------
      if (sensor % sens_type == SC_GMM_ID) then
         if (controlVariables % containsKey(SC_NUM_CLUSTERS_KEY)) then
            sensor % nclusters = controlVariables % doublePrecisionValueForKey(SC_NUM_CLUSTERS_KEY)
         else
            sensor % nclusters = 2
         end if
!
!     Options for the rest of sensors
!     -------------------------------
      else
!
!         Sensor thresholds
!         --------------
          if (controlVariables % containsKey(SC_LOW_THRES_KEY)) then
             sensor % low = controlVariables % doublePrecisionValueForKey(SC_LOW_THRES_KEY)
          else
             write(STD_OUT,*) 'ERROR. Lower threshold of the sensor must be specified.'
             stop
          end if

          if (controlVariables % containsKey(SC_HIGH_THRES_KEY)) then
             sensor % high = controlVariables % doublePrecisionValueForKey(SC_HIGH_THRES_KEY)
          else
             write(STD_OUT,*) 'ERROR. Higher threshold of the sensor must be specified.'
             stop
          end if
!
!         Sensor parameters
!         -----------------
          sensor % s0   = (sensor % high + sensor % low) / 2.0_RP
          sensor % ds   = (sensor % high - sensor % low)
          sensor % ds2  = sensor % ds / 2.0_RP

      end if
!
!     Sensed variable
!     ---------------
      if (sensor % sens_type == SC_GMM_ID) then
         sensor % sVar = SC_RHO_P_GRAD_ID

      elseif (controlVariables % containsKey(SC_VARIABLE_KEY)) then
         sVar = controlVariables % stringValueForKey(SC_VARIABLE_KEY, LINE_LENGTH)
         call toLower(sVar)

         select case (trim(sVar))
         case (SC_RHO_VAL);      sensor % sVar = SC_RHO_ID
         case (SC_RHOU_VAL);     sensor % sVar = SC_RHOU_ID
         case (SC_RHOV_VAL);     sensor % sVar = SC_RHOV_ID
         case (SC_RHOW_VAL);     sensor % sVar = SC_RHOW_ID
         case (SC_RHOE_VAL);     sensor % sVar = SC_RHOE_ID
         case (SC_U_VAL);        sensor % sVar = SC_U_ID
         case (SC_V_VAL);        sensor % sVar = SC_V_ID
         case (SC_W_VAL);        sensor % sVar = SC_W_ID
         case (SC_P_VAL);        sensor % sVar = SC_P_ID
         case (SC_RHOP_VAL);     sensor % sVar = SC_RHOP_ID
         case (SC_RHO_GRAD_VAL); sensor % sVar = SC_RHO_GRAD_ID
         case default
            write(STD_OUT,*) 'ERROR. The sensor variable is unknown. Options are:'
            write(STD_OUT,*) '   * ', SC_RHO_VAL
            write(STD_OUT,*) '   * ', SC_RHOU_VAL
            write(STD_OUT,*) '   * ', SC_RHOV_VAL
            write(STD_OUT,*) '   * ', SC_RHOW_VAL
            write(STD_OUT,*) '   * ', SC_RHOE_VAL
            write(STD_OUT,*) '   * ', SC_U_VAL
            write(STD_OUT,*) '   * ', SC_V_VAL
            write(STD_OUT,*) '   * ', SC_W_VAL
            write(STD_OUT,*) '   * ', SC_P_VAL
            write(STD_OUT,*) '   * ', SC_RHOP_VAL
            write(STD_OUT,*) '   * ', SC_RHO_GRAD_VAL
            errorMessage(STD_OUT)
            stop
         end select

      else
         sensor % sVar = SC_RHOP_ID

      end if

      if (sensor % sVar == SC_RHO_GRAD_ID) then
         if (.not. grad_vars == GRADVARS_STATE .and. .not. grad_vars == GRADVARS_ENTROPY) then
            write(STD_OUT,*) "ERROR. The density gradient sensor only works with state or entropy variables."
            stop
         end if
      end if

   end subroutine Set_SCsensor
!
!///////////////////////////////////////////////////////////////////////////////
!
   subroutine Construct_TEsensor(sensor, controlVariables, sem, &
                                 ComputeTimeDerivative, ComputeTimeDerivativeIsolated)
!
!     -------
!     Modules
!     -------
      use FTValueDictionaryClass
      use TransfiniteMapClass,        only: TransfiniteHexMap
      use NodalStorageClass,          only: NodalStorage
      use SpectralVanishingViscosity, only: SVV
!
!     ---------
!     Interface
!     ---------
      implicit none
      type(SCsensor_t),                   intent(inout) :: sensor
      type(FTValueDictionary),            intent(in)    :: controlVariables
      type(DGSem),                        intent(in)    :: sem
      procedure(ComputeTimeDerivative_f)                :: ComputeTimeDerivative
      procedure(ComputeTimeDerivative_f)                :: ComputeTimeDerivativeIsolated
!
!     ---------------
!     Local variables
!     ---------------
      character(len=:), allocatable    :: derivType
      real(RP)                         :: Nmin
      real(RP)                         :: deltaN
      integer                          :: eID
      integer,          allocatable    :: N(:,:)

!
!     Read the control file
!     ---------------------
      allocate(sensor % TEestim)
      if (controlVariables % containsKey(SC_TE_NMIN_KEY)) then
         Nmin = controlVariables % integerValueForKey(SC_TE_NMIN_KEY)
      else
         Nmin = 1
      end if

      if (controlVariables % containsKey(SC_TE_DELTA_KEY)) then
         deltaN = controlVariables % integerValueForKey(SC_TE_DELTA_KEY)
      else
         deltaN = 1
      end if

      if (controlVariables % containsKey(SC_TE_DTYPE_KEY)) then
         derivType = controlVariables % stringValueForKey(SC_TE_DTYPE_KEY, LINE_LENGTH)
      else
         derivType = SC_ISOLATED_KEY
      end if
      call toLower(derivType)

      select case (derivType)
      case (SC_ISOLATED_KEY)
         sensor % TEestim % derivType = ISOLATED_TE
         sensor % TEestim % TimeDerivative => ComputeTimeDerivativeIsolated
      case (SC_NON_ISOLATED_KEY)
         sensor % TEestim % derivType = NON_ISOLATED_TE
         sensor % TEestim % TimeDerivative => ComputeTimeDerivative
      case default
         write(STD_OUT,*) "ERROR. The TE sensor can only use the (non-)isolated time derivatives."
         stop
      end select

      sensor % TEestim % Nmin   = Nmin
      sensor % TEestim % deltaN = deltaN
!
!     Coarse mesh
!     -----------
      allocate(sensor % TEestim % coarseSem)
      sensor % TEestim % coarseSem = sem
      sensor % TEestim % coarseSem % mesh % child = .true.
      sensor % TEestim % coarseSem % mesh % ignoreBCnonConformities = .true.

      allocate(N(NDIM, sem % mesh % no_of_elements))
!$omp parallel do schedule(runtime)
      do eID = 1, sem % mesh % no_of_elements

         if (sem % mesh % Nx(eID) < Nmin) then
            N(IX,eID) = sem % mesh % Nx(eID)
         else
            N(IX,eID) = max(sem % mesh % Nx(eID) - deltaN, Nmin)
         end if

         if (sem % mesh % Ny(eID) < Nmin) then
            N(IY,eID) = sem % mesh % Ny(eID)
         else
            N(IY,eID) = max(sem % mesh % Ny(eID) - deltaN, Nmin)
         end if

         if (sem % mesh % Nz(eID) < Nmin) then
            N(IZ,eID) = sem % mesh % Nz(eID)
         else
            N(IZ,eID) = max(sem % mesh % Nz(eID) - deltaN, Nmin)
         end if

      end do
!$omp end parallel do

      call sensor % TEestim % coarseSem % mesh % pAdapt(N, controlVariables)
      call sensor % TEestim % coarseSem % mesh % storage % PointStorage()
!
!     Update the SVV if active
!     ------------------------
      if (SVV % enabled) then
         call SVV % UpdateFilters(sensor % TEestim % coarseSem % mesh)
      end if

   end subroutine Construct_TEsensor
!
!///////////////////////////////////////////////////////////////////////////////
!
   subroutine Describe_SCsensor(sensor)
!
!     ---------
!     Interface
!     ---------
      implicit none
      class(SCsensor_t), intent(in) :: sensor


      write(STD_OUT,"(30X,A,A30)", advance="no") "->", "Sensor type: "
      select case (sensor % sens_type)
         case (SC_ZERO_ID);     write(STD_OUT,"(A)") SC_ZERO_VAL
         case (SC_ONE_ID);      write(STD_OUT,"(A)") SC_ONE_VAL
         case (SC_INTEGRAL_ID); write(STD_OUT,"(A)") SC_INTEGRAL_VAL
         case (SC_MODAL_ID);    write(STD_OUT,"(A)") SC_MODAL_VAL
         case (SC_TE_ID);       write(STD_OUT,"(A)") SC_TE_VAL
         case (SC_ALIAS_ID);    write(STD_OUT,"(A)") SC_ALIAS_VAL
         case (SC_GMM_ID);      write(STD_OUT,"(A)") SC_GMM_VAL
      end select

      if (sensor % sens_type == SC_TE_ID) then
         write(STD_OUT,"(30X,A,A30)", advance="no") "->", "Time derivative: "
         select case (sensor % TEestim % derivType)
            case (ISOLATED_TE);     write(STD_OUT,"(A)") SC_ISOLATED_KEY
            case (NON_ISOLATED_TE); write(STD_OUT,"(A)") SC_NON_ISOLATED_KEY
         end select

         write(STD_OUT,"(30X,A,A30,I0,A,I0,A)") "->", "Delta N: ", sensor % TEestim % deltaN, &
                                                " (min. order is ", sensor % TEestim % Nmin, ")"
      end if

      write(STD_OUT,"(30X,A,A30)", advance="no") "->", "Sensed variable: "
      select case (sensor % sVar)
         case (SC_RHO_ID);        write(STD_OUT,"(A)") SC_RHO_VAL
         case (SC_RHOU_ID);       write(STD_OUT,"(A)") SC_RHOU_VAL
         case (SC_RHOV_ID);       write(STD_OUT,"(A)") SC_RHOV_VAL
         case (SC_RHOW_ID);       write(STD_OUT,"(A)") SC_RHOW_VAL
         case (SC_RHOE_ID);       write(STD_OUT,"(A)") SC_RHOE_VAL
         case (SC_U_ID);          write(STD_OUT,"(A)") SC_U_VAL
         case (SC_V_ID);          write(STD_OUT,"(A)") SC_V_VAL
         case (SC_W_ID);          write(STD_OUT,"(A)") SC_W_VAL
         case (SC_P_ID);          write(STD_OUT,"(A)") SC_P_VAL
         case (SC_RHOP_ID);       write(STD_OUT,"(A)") SC_RHOP_VAL
         case (SC_RHO_GRAD_ID);   write(STD_OUT,"(A)") SC_RHO_GRAD_VAL
         case (SC_RHO_P_GRAD_ID); write(STD_OUT,"(A)") SC_RHO_P_GRAD_VAL
      end select

      write(STD_OUT,"(30X,A,A30,I0,A)") "->", "Sensor inertia: ", sensor % min_steps, " timesteps"

      if (sensor % sens_type == SC_GMM_ID) then
         write(STD_OUT,"(30X,A,A30,I0)") "->", "Number of clusters: ", sensor % nclusters
      else
         write(STD_OUT,"(30X,A,A30,F5.1)") "->", "Minimum value: ", sensor % low
         write(STD_OUT,"(30X,A,A30,F5.1)") "->", "Maximum value: ", sensor % high
      end if


   end subroutine Describe_SCsensor
!
!///////////////////////////////////////////////////////////////////////////////
!
   subroutine Destruct_SCsensor(sensor)
!
!     ---------
!     Interface
!     ---------
      type(SCsensor_t), intent(inout) :: sensor


      if (allocated(sensor % TEestim))      deallocate(sensor % TEestim)
      if (associated(sensor % Compute_Raw)) nullify(sensor % Compute_Raw)

   end subroutine Destruct_SCsensor
!
!///////////////////////////////////////////////////////////////////////////////
!
   subroutine Compute_SCsensor(sensor, sem, t)
!
!     ---------
!     Interface
!     ---------
      implicit none
      class(SCsensor_t), target, intent(inout) :: sensor
      type(DGsem),       target, intent(inout) :: sem
      real(RP),                  intent(in)    :: t
!
!     ---------------
!     Local variables
!     ---------------
!
      integer  :: eID

!
!     Compute the sensor
!     ------------------
      call sensor % Compute_Raw(sem, t)
!
!     Add 'inertia' to the scaled value
!     ---------------------------------
!$omp parallel do
      do eID = 1, sem % mesh % no_of_elements
      associate(e     => sem % mesh % elements(eID),                  &
                value => sem % mesh % elements(eID) % storage % sensor)
         if (value <= 0.0_RP .and. e % storage % last_sensed < sensor % min_steps) then
            e % storage % sensor = e % storage % prev_sensor
            e % storage % last_sensed = e % storage % last_sensed + 1
         else
            e % storage % sensor = value
            e % storage % prev_sensor = value
            e % storage % last_sensed = 0
         end if
      end associate
      end do
!$omp end parallel do

   end subroutine Compute_SCsensor
!
!///////////////////////////////////////////////////////////////////////////////
!
!     Sensors
!
!///////////////////////////////////////////////////////////////////////////////
!
   subroutine Sensor_zero(sensor, sem, t)
!
!     ---------
!     Interface
!     ---------
      implicit none
      class(SCsensor_t), target, intent(inout) :: sensor
      type(DGsem),       target, intent(inout) :: sem
      real(RP),                  intent(in)    :: t
!
!     ---------------
!     Local variables
!     ---------------
      integer :: eID


!$omp parallel do
      do eID = 1, sem % mesh % no_of_elements
         sem % mesh % elements(eID) % storage % sensor = 0.0_RP
      end do
!$omp end parallel do

   end subroutine Sensor_zero
!
!///////////////////////////////////////////////////////////////////////////////
!
   subroutine Sensor_one(sensor, sem, t)
!
!     ---------
!     Interface
!     ---------
      implicit none
      class(SCsensor_t), target, intent(inout) :: sensor
      type(DGSem),       target, intent(inout) :: sem
      real(RP),                  intent(in)    :: t
!
!     ---------------
!     Local variables
!     ---------------
      integer :: eID


!$omp parallel do
      do eID = 1, sem % mesh % no_of_elements
         sem % mesh % elements(eID) % storage % sensor = 1.0_RP
      end do
!$omp end parallel do

   end subroutine Sensor_one
!
!///////////////////////////////////////////////////////////////////////////////
!
   subroutine Sensor_integral(sensor, sem, t)
!
!     -------
!     Modules
!     -------
      use NodalStorageClass, only: NodalStorage
!
!     ---------
!     Interface
!     ---------
      implicit none
      class(SCsensor_t), target, intent(inout) :: sensor
      type(DGSem),       target, intent(inout) :: sem
      real(RP),                  intent(in)    :: t
!
!     ---------------
!     Local variables
!     ---------------
      integer  :: eID
      integer  :: i
      integer  :: j
      integer  :: k
      real(RP) :: contribution
      real(RP) :: val

!$omp parallel do schedule(runtime)
      do eID = 1, sem % mesh % no_of_elements
      associate(e => sem % mesh % elements(eID))

         val = 0.0_RP
         do k = 0, e % Nxyz(3) ; do j = 0, e % Nxyz(2) ; do i = 0, e % Nxyz(1)
            contribution = GetSensedVariable( &
               sensor % sVar,                 &
               e % storage % Q(:,i,j,k),      &
               e % storage % U_x(:,i,j,k),    &
               e % storage % U_y(:,i,j,k),    &
               e % storage % U_z(:,i,j,k)     &
            )
            val = val + NodalStorage(e % Nxyz(1)) % w(i) &
                      * NodalStorage(e % Nxyz(2)) % w(j) &
                      * NodalStorage(e % Nxyz(3)) % w(k) &
                      * e % geom % jacobian(i,j,k)       &
                      * contribution
         end do                ; end do                ; end do

         e % storage % sensor = SinRamp(sensor, sqrt(val))

      end associate
      end do
!$omp end parallel do

   end subroutine Sensor_integral
!
!///////////////////////////////////////////////////////////////////////////////
!
   subroutine Sensor_modal(sensor, sem, t)
!
!     -------
!     Modules
!     -------
      use Utilities, only: AlmostEqual
!
!     ---------
!     Interface
!     ---------
      implicit none
      class(SCsensor_t), target, intent(inout) :: sensor
      type(DGSem),       target, intent(inout) :: sem
      real(RP),                  intent(in)    :: t
!
!     ---------------
!     Local variables
!     ---------------
      integer               :: eID
      integer               :: i, j, k, r
      integer               :: maxNx, maxNy, maxNz
      real(RP), allocatable :: sVar(:,:,:)
      real(RP), allocatable :: sVarMod(:,:,:)
      real(RP), allocatable :: Lwx(:), Lwy(:), Lwz(:)
      real(RP)              :: num, den


      ! Only allocate once
      maxNx = maxval(sem % mesh % Nx)
      maxNy = maxval(sem % mesh % Ny)
      maxNz = maxval(sem % mesh % Nz)
      allocate(sVar(maxNx,maxNy,maxNz))
      allocate(sVarMod(maxNx,maxNy,maxNz))

!$omp parallel do schedule(runtime) private(eID, sVar, sVarMod)
      do eID = 1, sem % mesh % no_of_elements
      associate(e => sem % mesh % elements(eID))

         associate(Nx => e % Nxyz(1), Ny => e % Nxyz(2), Nz => e % Nxyz(3))
!
!        Compute the sensed variable
!        ---------------------------
         do k = 0, Nz ; do j = 0, Ny ; do i = 0, Nx
            sVar(i,j,k) = GetSensedVariable( &
               sensor % sVar,                &
               e % storage % Q(:,i,j,k),     &
               e % storage % U_x(:,i,j,k),   &
               e % storage % U_y(:,i,j,k),   &
               e % storage % U_z(:,i,j,k)    &
            )
         end do       ; end do       ; end do
!
!        Switch to modal space
!        ---------------------
         sVarMod = 0.0_RP

         ! Xi direction
         do k = 0, Nz ; do j = 0, Ny ; do r = 0, Nx ; do i = 0, Nx
            sVarMod(i,j,k) = sVarMod(i,j,k) + NodalStorage(Nx) % Fwd(i,r) * sVar(r,j,k)
         end do       ; end do       ; end do       ; end do

         ! Eta direction
         do k = 0, Nz ; do r = 0, Ny ; do j = 0, Ny ; do i = 0, Nx
            sVarMod(i,j,k) = sVarMod(i,j,k) + NodalStorage(Ny) % Fwd(j,r) * sVar(i,r,k)
         end do       ; end do       ; end do       ; end do

         ! Zeta direction
         do r = 0, Nz ; do k = 0, Nz ; do j = 0, Ny ; do i = 0, Nx
            sVarMod(i,j,k) = sVarMod(i,j,k) + NodalStorage(Nz) % Fwd(k,r) * sVar(i,j,r)
         end do       ; end do       ; end do       ; end do
!
!        Check almost zero values
!        ------------------------
         do k = 0, Nz ; do j = 0, Ny ; do i = 0 , Nx
            if (AlmostEqual(sVarMod(i,j,k), 0.0_RP)) then
               sVarMod(i,j,k) = abs(sVarMod(i,j,k))
            end if
         end do       ; end do       ; end do
!
!        Ratio of higher modes vs all the modes
!        --------------------------------------
!        Explanation: The higher modes are contained in the subspace V(Nx,Ny,Nz) - V(Nx-1,Ny-1,Nz-1)
!                     Representing them as cubes, this subspace corresponds to the outer shell of
!                     V(Nx,Ny,Nz), thus the terms "faces", "edges" and "corner"
!        --------------------------------------------------------------------------------------------
         Lwx = NodalStorage(Nx) % Lw
         Lwy = NodalStorage(Ny) % Lw
         Lwz = NodalStorage(Nz) % Lw

         ! Corner (Nx, Ny, Nz)
         num = sVarMod(Nx,Ny,Nz)**2 * Lwx(Nx) * Lwy(Ny) * Lwz(Nz)

         ! Edges
         do i = 0, Nx-1  ! +X edge
            num = num + sVarMod(i,Ny,Nz)**2 * Lwx(i) * Lwy(Ny) * Lwz(Nz)
         end do
         do j = 0, Ny-1  ! +Y edge
            num = num + sVarMod(Nx,j,Nz)**2 * Lwx(Nx) * Lwy(j) * Lwz(Nz)
         end do
         do k = 0, Nz-1  ! +Z edge
            num = num + sVarMod(Nx,Ny,k)**2 * Lwx(Nx) * Lwy(Ny) * Lwz(k)
         end do

         ! Faces
         do k = 0, Nz-1 ; do j = 0, Ny-1  ! +X face
            num = num + sVarMod(Nx,j,k)**2 * Lwx(Nx) * Lwy(j) * Lwz(k)
         end do         ; end do
         do k = 0, Nz-1 ; do i = 0, Nx-1  ! +Y face
            num = num + sVarMod(i,Ny,k)**2 * Lwx(i) * Lwy(Ny) * Lwz(k)
         end do         ; end do
         do j = 0, Ny-1 ; do i = 0, Nx-1  ! +Z face
            num = num + sVarMod(i,j,Nz)**2 * Lwx(i) * Lwy(j) * Lwz(Nz)
         end do         ; end do

         ! Total sum
         den = 0.0_RP
         do k = 0, Nz ; do j = 0, Ny ; do i = 0, Nx
            den = den + sVarMod(i,j,k)**2 * Lwx(i) * Lwy(j) * Lwz(k)
         end do       ; end do       ; end do

         ! Sensor value as the ratio of num / den
         e % storage % sensor = SinRamp(sensor, log10( num / den ))

         end associate

      end associate
      end do
!$omp end parallel do

   end subroutine Sensor_modal
!
!///////////////////////////////////////////////////////////////////////////////
!
   subroutine Sensor_truncation(sensor, sem, t)
!
!     -------
!     Modules
!     -------
      use InterpolationMatrices, only: Interp3DArrays
      use ProblemFileFunctions,  only: UserDefinedSourceTermNS_f
      use PhysicsStorage,        only: CTD_IGNORE_MODE
      use FluidData,             only: thermodynamics, dimensionless, refValues
      use Utilities,             only: AlmostEqual
      use ElementClass,          only: Element
!
!     ---------
!     Interface
!     ---------
      implicit none
      class(SCsensor_t), target, intent(inout) :: sensor
      type(DGSem),       target, intent(inout) :: sem
      real(RP),                  intent(in)    :: t
!
!     ---------------
!     Local variables
!     ---------------
      integer                              :: eID
      integer                              :: i, j, k
      integer                              :: iEq
      real(RP)                             :: wx, wy, wz
      real(RP)                             :: Jac
      real(RP)                             :: S(NCONS)
      real(RP)                             :: mTE
      procedure(UserDefinedSourceTermNS_f) :: UserDefinedSourceTermNS

      ! TODO: ifort does not like associates with OpenMP
      type(DGSem),   pointer :: csem => null()
      type(Element), pointer :: e  => null()
      type(Element), pointer :: ce => null()


      csem => sensor % TEestim % coarseSem
!
!     Time derivative
!     ---------------
!$omp parallel do schedule(runtime) private(eID, e, ce)
      do eID = 1, sem % mesh % no_of_elements

         e  => sem % mesh % elements(eID)
         ce => csem % mesh % elements(eID)

         call Interp3DArrays(NCONS, e % Nxyz, e % storage % Q, ce % Nxyz, ce % storage % Q)
         ce % storage % sensor = e % storage % sensor

      end do
!$omp end parallel do

      call sensor % TEestim % TimeDerivative(csem % mesh, csem % particles, t, CTD_IGNORE_MODE)
!
!     Maximum TE computation
!     ----------------------
!$omp parallel do schedule(runtime) private(eID, e, ce, i, j, k, iEq, S, mTE, wx, wy, wz, Jac)
      do eID = 1, sem % mesh % no_of_elements

         e  => sem % mesh % elements(eID)
         ce => csem % mesh % elements(eID)

         ! Use G_NS as temporary storage for Qdot
         call Interp3DArrays(NCONS, e % Nxyz, e % storage % QDot, ce % Nxyz, ce % storage % G_NS)

         mTE = 0.0_RP
         do k = 0, ce % Nxyz(IZ) ; do j = 0, ce % Nxyz(IY) ; do i = 0, ce % Nxyz(IX)

            call UserDefinedSourceTermNS(ce % geom % x, ce % storage % Q(:,i,j,k), t, S, &
                                         thermodynamics, dimensionless, refValues)
            wx = NodalStorage(ce % Nxyz(IX)) % w(i)
            wy = NodalStorage(ce % Nxyz(IY)) % w(j)
            wz = NodalStorage(ce % Nxyz(IZ)) % w(k)
            Jac = ce % geom % jacobian(i,j,k)

            mTE = max(mTE, Jac*wx*wy*wz * maxval(abs(ce % storage % QDot(:,i,j,k) + S(:) - &
                                                     ce % storage % G_NS(:,i,j,k))))


         end do                  ; end do                  ; end do

         if (AlmostEqual(mTE, 0.0_RP)) then
            e % storage % sensor = 0.0_RP
         else
            e % storage % sensor = SinRamp(sensor, log10(mTE))
         end if

      end do
!$omp end parallel do

      nullify(csem)
      if (associated(e))  nullify(e)
      if (associated(ce)) nullify(ce)

   end subroutine Sensor_truncation
!
!///////////////////////////////////////////////////////////////////////////////
!
   subroutine Sensor_aliasing(sensor, sem, t)
!
!     -------
!     Modules
!     -------
      use HyperbolicDiscretizations, only: HyperbolicDiscretization
      use Physics,                   only: EulerFlux
      use Utilities,                 only: AlmostEqual
!
!     ---------
!     Interface
!     ---------
      implicit none
      class(SCsensor_t), target, intent(inout) :: sensor
      type(DGSem),       target, intent(inout) :: sem
      real(RP),                  intent(in)    :: t
!
!     ---------------
!     Local variables
!     ---------------
      integer               :: eID
      integer               :: i, j, k, l
      real(RP), allocatable :: F(:,:,:,:,:)
      real(RP), allocatable :: Fs(:,:,:,:,:)
      real(RP), allocatable :: Gs(:,:,:,:,:)
      real(RP), allocatable :: Hs(:,:,:,:,:)
      real(RP), allocatable :: aliasing(:,:,:,:)


!$omp parallel do schedule(runtime) private(eID, i, j, k, l, F, Fs, Gs, Hs, aliasing)
      do eID = 1, sem % mesh % no_of_elements
      associate(e => sem % mesh % elements(eID))

         allocate(F(NCONS,  0:e % Nxyz(1), 0:e % Nxyz(2), 0:e % Nxyz(3), NDIM))
         allocate(Fs(NCONS, 0:e % Nxyz(1), 0:e % Nxyz(1), 0:e % Nxyz(2), 0:e % Nxyz(3)))
         allocate(Gs(NCONS, 0:e % Nxyz(2), 0:e % Nxyz(1), 0:e % Nxyz(2), 0:e % Nxyz(3)))
         allocate(Hs(NCONS, 0:e % Nxyz(3), 0:e % Nxyz(1), 0:e % Nxyz(2), 0:e % Nxyz(3)))

         call HyperbolicDiscretization % ComputeInnerFluxes(e, EulerFlux, F)
         call HyperbolicDiscretization % ComputeSplitFormFluxes(e, F, Fs, Gs, Hs)

         associate(spAxi   => NodalStorage(e % Nxyz(1)), &
                   spAeta  => NodalStorage(e % Nxyz(2)), &
                   spAzeta => NodalStorage(e % Nxyz(3)))

            allocate(aliasing(NCONS, 0:e % Nxyz(1), 0:e % Nxyz(2), 0:e % Nxyz(3)), source=0.0_RP)
            do k = 0, e % Nxyz(3) ; do j = 0, e % Nxyz(2) ; do l = 0, e % Nxyz(1) ; do i = 0, e % Nxyz(1)
               aliasing(:,i,j,k) = aliasing(:,i,j,k) + spAxi % D(i,l) * (2.0_RP*Fs(:,i,l,j,k) - Fs(:,l,l,j,k))
            end do                ; end do                ; end do                ; end do
            do k = 0, e % Nxyz(3) ; do l = 0, e % Nxyz(2) ; do j = 0, e % Nxyz(2) ; do i = 0, e % Nxyz(1)
               aliasing(:,i,j,k) = aliasing(:,i,j,k) + spAeta % D(j,l) * (2.0_RP*Gs(:,j,i,l,k) - Gs(:,l,i,l,k))
            end do                ; end do                ; end do                ; end do
            do l = 0, e % Nxyz(3) ; do k = 0, e % Nxyz(3) ; do j = 0, e % Nxyz(2) ; do i = 0, e % Nxyz(1)
               aliasing(:,i,j,k) = aliasing(:,i,j,k) + spAzeta % D(k,l) * (2.0_RP*Hs(:,k,i,j,l) - Hs(:,l,i,j,l))
            end do                ; end do                ; end do                ; end do

         end associate

         aliasing = abs(aliasing)
         e % storage % sensor = maxval(aliasing)
         if (AlmostEqual(e % storage % sensor, 0.0_RP)) then
            e % storage % sensor = 0.0_RP
         else
            e % storage % sensor = SinRamp(sensor, log10(e % storage % sensor))
         end if

         deallocate(F)
         deallocate(Fs)
         deallocate(Gs)
         deallocate(Hs)
         deallocate(aliasing)

      end associate
      end do
!$omp end parallel do

   end subroutine Sensor_aliasing
!
!///////////////////////////////////////////////////////////////////////////////
!
   subroutine Sensor_GMM(sensor, sem, t)
!
!     -------
!     Modules
!     -------
      use PhysicsStorage,     only: IRHO, IRHOU, IRHOV, IRHOW, IRHOE
      use FluidData,          only: thermodynamics
      use VariableConversion, only: getVelocityGradients
!
!     ---------
!     Interface
!     ---------
      implicit none
      class(SCsensor_t), target, intent(inout) :: sensor
      type(DGSem),       target, intent(inout) :: sem
      real(RP),                  intent(in)    :: t
!
!     ---------------
!     Local variables
!     ---------------
      integer  :: eID
      integer  :: i, j, k
      integer  :: ndofs
      integer  :: nclusters
      integer  :: cnt
      integer  :: n
      integer  :: higherCluster
      real(RP) :: u2
      real(RP) :: ux(3), uy(3), uz(3)
      real(RP) :: dp(3)
      real(RP), allocatable :: centroids(:,:)
      real(RP), allocatable :: derivs(:,:)
      integer,  allocatable :: clusters(:)


      ndofs = sem % mesh % NDOF
      allocate(centroids(2,sensor % nclusters))
      allocate(derivs(2,ndofs))
      allocate(clusters(ndofs))
!
!     Compute the clustering variables and store them in a global array
!     -----------------------------------------------------------------
      cnt = 0
      do eID = 1, sem % mesh % no_of_elements
      associate(e => sem % mesh % elements(eID))

         do k = 0, e % Nxyz(3) ; do j = 0, e % Nxyz(2) ; do i = 0, e % Nxyz(1)
            cnt = cnt + 1
            u2 = (                               &
               e % storage % Q(IRHOU,i,j,k)**2 + &
               e % storage % Q(IRHOV,i,j,k)**2 + &
               e % storage % Q(IRHOW,i,j,k)**2   &
            ) / e % storage % Q(IRHO,i,j,k)**2

            ! Derivative of density
            derivs(1,cnt) = sqrt(                 &
               e % storage % U_x(IRHO,i,j,k)**2 + &
               e % storage % U_y(IRHO,i,j,k)**2 + &
               e % storage % U_z(IRHO,i,j,k)**2   &
            )

            ! Derivative of pressure
            call getVelocityGradients(     &
               e % storage % Q(:,i,j,k),   &
               e % storage % U_x(:,i,j,k), &
               e % storage % U_y(:,i,j,k), &
               e % storage % U_z(:,i,j,k), &
               ux, uy, uz                  &
            )
            dp(1) = thermodynamics % gammaMinus1 * (         &
               e % storage % U_x(IRHOE,i,j,k) -              &
               0.5_RP * e % storage % U_x(IRHO,i,j,k) * u2 - &
               e % storage % Q(IRHOU,i,j,k) * ux(1) -        &
               e % storage % Q(IRHOV,i,j,k) * ux(2) -        &
               e % storage % Q(IRHOW,i,j,k) * ux(3)          &
            )
            dp(2) = thermodynamics % gammaMinus1 * (         &
               e % storage % U_y(IRHOE,i,j,k) -              &
               0.5_RP * e % storage % U_y(IRHO,i,j,k) * u2 - &
               e % storage % Q(IRHOU,i,j,k) * uy(1) -        &
               e % storage % Q(IRHOV,i,j,k) * uy(2) -        &
               e % storage % Q(IRHOW,i,j,k) * uy(3)          &
            )
            dp(3) = thermodynamics % gammaMinus1 * (         &
               e % storage % U_z(IRHOE,i,j,k) -              &
               0.5_RP * e % storage % U_z(IRHO,i,j,k) * u2 - &
               e % storage % Q(IRHOU,i,j,k) * uz(1) -        &
               e % storage % Q(IRHOV,i,j,k) * uz(2) -        &
               e % storage % Q(IRHOW,i,j,k) * uz(3)          &
            )
            derivs(2,cnt) = sqrt(dp(1)**2 + dp(2)**2 + dp(3)**2)
         end do ;                end do ;                end do

      end associate
      end do
!
!     Rescale the values
!     ------------------
      call GetClusterVariables(sensor % nclusters, [1.0_RP, 20.0_RP], derivs, centroids)
!
!     Compute the GMM clusters
!     ------------------------
      call GMM(sensor % nclusters, nclusters, derivs, centroids, clusters)
!
!     Compute the sensor values
!     -------------------------
      cnt = 0
      do eID = 1, sem % mesh % no_of_elements
      associate(e => sem % mesh % elements(eID))
         if (nclusters <= 1) then
            e % storage % sensor = 0.0_RP
         else
            n = product(e % Nxyz + 1)
            higherCluster = maxval(clusters(cnt+1:cnt+n))
            e % storage % sensor = real(higherCluster - 1, RP) / (nclusters - 1)
         end if
         cnt = cnt + n
      end associate
      end do

   end subroutine Sensor_GMM
!
!///////////////////////////////////////////////////////////////////////////////
!
!     Scaling functions
!
!///////////////////////////////////////////////////////////////////////////////
!
   pure function SinRamp(sensor, sensedValue) result(scaled)
!
!     ---------
!     Interface
!     ---------
      implicit none
      class(SCsensor_t), intent(in) :: sensor
      real(RP),          intent(in) :: sensedValue
      real(RP)                      :: scaled
!
!     ---------------
!     Local variables
!     ---------------
      real(RP) :: dev


      ! Deviation of the input value wrt the centre of the sensor
      dev = sensedValue - sensor % s0

      ! Sin-ramp activation function
      if (dev <= -sensor % ds2) then
         scaled = 0.0_RP

      else if (dev >= sensor % ds2) then
         scaled = 1.0_RP

      else
         scaled = ( 1.0_RP + sin(PI * dev / sensor % ds) ) / 2.0_RP

      end if

   end function SinRamp
!
!///////////////////////////////////////////////////////////////////////////////
!
!     Utilities
!
!///////////////////////////////////////////////////////////////////////////////
!
   pure function GetSensedVariable(varType, Q, U_x, U_y, U_z) result(s)
!
!     -------
!     Modules
!     -------
      use PhysicsStorage,     only: IRHO, IRHOU, IRHOV, IRHOW, IRHOE
      use VariableConversion, only: Pressure
      implicit none
!
!     ---------
!     Interface
!     ---------
      integer,  intent(in) :: varType
      real(RP), intent(in) :: Q(:)
      real(RP), intent(in) :: U_x(:)
      real(RP), intent(in) :: U_y(:)
      real(RP), intent(in) :: U_z(:)
      real(RP)             :: s


      select case (varType)
      case (SC_RHO_ID)
         s = Q(IRHO)

      case (SC_RHOU_ID)
         s = Q(IRHOU)

      case (SC_RHOV_ID)
         s = Q(IRHOV)

      case (SC_RHOW_ID)
         s = Q(IRHOW)

      case (SC_RHOE_ID)
         s = Q(IRHOE)

      case (SC_U_ID)
         s = Q(IRHOU) / Q(IRHO)

      case (SC_V_ID)
         s = Q(IRHOV) / Q(IRHO)

      case (SC_W_ID)
         s = Q(IRHOW) / Q(IRHO)

      case (SC_P_ID)
         s = Pressure(Q)

      case (SC_RHOP_ID)
         s = Pressure(Q) * Q(IRHO)

      case (SC_RHO_GRAD_ID)
         if ( grad_vars == GRADVARS_STATE ) then
            s = POW2(U_x(IRHO)) + POW2(U_y(IRHO)) + POW2(U_z(IRHO))
         else if (grad_vars == GRADVARS_ENTROPY) then
            s = POW2(sum(Q * U_x)) + POW2(sum(Q * U_y)) + POW2(sum(Q * U_z))
         end if

      end select

   end function GetSensedVariable
!
!///////////////////////////////////////////////////////////////////////////////
!
   pure subroutine GetClusterVariables(nclusters, limit, x, xavg)
!
!     -------
!     Modules
!     -------
      use Utilities, only: AlmostEqual
!
!     ---------
!     Interface
!     ---------
      integer,  intent(in)    :: nclusters
      real(RP), intent(in)    :: limit(:)
      real(RP), intent(inout) :: x(:,:)
      real(RP), intent(out)   :: xavg(:,:)
!
!     ---------------
!     Local variables
!     ---------------
      integer  :: i
      real(RP) :: minimum(2)
      real(RP) :: maximum(2)
      real(RP) :: diff(2)


      x = abs(x)
      minimum = minval(x, dim=2)
      maximum = maxval(x, dim=2)

      if (AlmostEqual(maximum(1), minimum(1))) then
         if (maximum(1) > 0.0_RP) then
            x(1,:) = limit(1)
            minimum(1) = limit(1)
            maximum(1) = limit(1)
         else
            x(1,:) = 0.0_RP
            minimum(1) = 0.0_RP
            maximum(1) = 0.0_RP
         end if
      else
         x(1,:) = (x(1,:) - minimum(1)) / (maximum(1) - minimum(1))
         minimum(1) = 0.0_RP
         maximum(1) = limit(1)
      end if

      if (AlmostEqual(maximum(2), minimum(2))) then
         if (maximum(2) > 0.0_RP) then
            x(2,:) = limit(2)
            minimum(2) = limit(2)
            maximum(2) = limit(2)
         else
            x(2,:) = 0.0_RP
            minimum(2) = 0.0_RP
            maximum(2) = 0.0_RP
         end if
      else
         x(2,:) = (x(2,:) - minimum(2)) / (maximum(2) - minimum(2))
         minimum(2) = 0.0_RP
         maximum(2) = limit(2)
      end if

      xavg(:,1) = minimum
      xavg(:,nclusters) = maximum
      diff = (maximum - minimum) / (nclusters - 1)
      do i = 2, nclusters - 1
         xavg(:,i) = minimum + (i - 1) * diff
      end do

   end subroutine GetClusterVariables

end module SCsensorClass
