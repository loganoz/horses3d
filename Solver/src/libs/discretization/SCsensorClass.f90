!
!//////////////////////////////////////////////////////
!
!   @File:    ShockCapturing.f90
!   @Author:  Andrés Mateo (andres.mgabin@upm.es)
!   @Created: Thu Jun 17 2021
!   @Last revision date: Thu Wed 28 2021
!   @Last revision author: Andrés Mateo (andres.mgabin@upm.es)
!
!//////////////////////////////////////////////////////
!
#include "Includes.h"
module SCsensorClass

   use SMConstants,       only: RP, PI, STD_OUT, NDIM, IX, IY, IZ, LINE_LENGTH
   use DGSEMClass,        only: ComputeTimeDerivative_f, DGSem
   use PhysicsStorage,    only: grad_vars, GRADVARS_STATE, GRADVARS_ENTROPY, NCONS
   use NodalStorageClass, only: NodalStorage
   use Utilities,         only: toLower

   use ShockCapturingKeywords

   implicit none

   private

   public :: SCsensor_t
   public :: Set_SCsensor
   public :: Destruct_SCsensor

   type :: SCsensor_t

         integer  :: sens_type  ! ID of the sensor type

         real(RP) :: s0    !< Centre of the sensor scaling
         real(RP) :: ds    !< Bandwith of the sensor scaling
         real(RP) :: ds2   !< Half-bandwidth
         real(RP) :: low   !< Lower threshold
         real(RP) :: high  !< Upper threshold
         integer  :: sVar  !< Variable used as input for the sensor

         type(TruncationError_t), allocatable :: TEestim  !< Truncation error estimation

         procedure(Compute_Int), pointer :: Compute => null()
         procedure(Rescale_Int), pointer :: Rescale => SinRamp

      contains

         procedure :: Describe => Describe_SCsensor
         final     :: Destruct_SCsensor

   end type SCsensor_t

   type :: TruncationError_t
      integer                  :: Nmin       !< Min order N
      integer                  :: deltaN     !< N' of coarser mesh is max(Nmin, N-deltaN)
      integer                  :: derivType  !< Whether the time derivative isolated or not
      type(DGSem), allocatable :: coarseSem  !< DGSEM for the coarser mesh
      procedure(ComputeTimeDerivative_f), nopass, pointer :: TimeDerivative => null()
   end type TruncationError_t
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

      pure function Rescale_Int(this, sensVal) result(scaled)
         import RP, SCsensor_t
         class(SCsensor_t), intent(in) :: this
         real(RP),          intent(in) :: sensVal
         real(RP)                      :: scaled
      end function Rescale_Int
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
   subroutine Set_SCsensor(sensor, controlVariables, sem, complementaryGeometry, &
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
      logical,                            intent(in)    :: complementaryGeometry
      procedure(ComputeTimeDerivative_f)                :: ComputeTimeDerivative
      procedure(ComputeTimeDerivative_f)                :: ComputeTimeDerivativeIsolated
!
!     ---------------
!     Local variables
!     ---------------
      character(len=:), allocatable :: sensorType
      character(len=:), allocatable :: sVar

!
!     Sensor type
!     -----------
      if (controlVariables % containsKey(SC_SENSOR_KEY)) then
         sensorType = controlVariables % stringValueForKey(SC_SENSOR_KEY, LINE_LENGTH)
      else
         sensorType = SC_GRADRHO_VAL
      end if
      call toLower(sensorType)

      select case (sensorType)
      case (SC_ZERO_VAL)
         sensor % sens_type = SC_ZERO_ID
         sensor % Compute  => Sensor_zero

      case (SC_ONE_VAL)
         sensor % sens_type = SC_ONE_ID
         sensor % Compute  => Sensor_one

      case (SC_GRADRHO_VAL)
         if (grad_vars == GRADVARS_STATE .or. grad_vars == GRADVARS_ENTROPY) then
            sensor % sens_type = SC_GRADRHO_ID
            sensor % Compute  => Sensor_rho
         else
            write(STD_OUT,*) "ERROR. The density sensor only works with state or entropy variables."
            stop
         end if

      case (SC_MODAL_VAL)
         sensor % sens_type = SC_MODAL_ID
         sensor % Compute  => Sensor_modal

      case (SC_TE_VAL)
         sensor % sens_type = SC_TE_ID
         sensor % Compute  => Sensor_truncation
         call Construct_TEsensor(sensor, controlVariables, sem, complementaryGeometry, &
                                 ComputeTimeDerivative, ComputeTimeDerivativeIsolated)

      case default
         write(STD_OUT,*) "ERROR. The sensor type is unkown. Options are:"
         write(STD_OUT,*) '   * ', SC_ZERO_VAL
         write(STD_OUT,*) '   * ', SC_ONE_VAL
         write(STD_OUT,*) '   * ', SC_GRADRHO_VAL
         write(STD_OUT,*) '   * ', SC_MODAL_VAL
         write(STD_OUT,*) '   * ', SC_TE_VAL
         stop

      end select
!
!     Sensor thresholds
!     --------------
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
!     Sensed variable
!     ---------------
      if (controlVariables % containsKey(SC_VARIABLE_KEY)) then
         sVar = controlVariables % stringValueForKey(SC_VARIABLE_KEY, LINE_LENGTH)
         call toLower(sVar)

         select case (trim(sVar))
         case (SC_RHO_VAL);  sensor % sVar = SC_RHO_ID
         case (SC_RHOU_VAL); sensor % sVar = SC_RHOU_ID
         case (SC_RHOV_VAL); sensor % sVar = SC_RHOV_ID
         case (SC_RHOW_VAL); sensor % sVar = SC_RHOW_ID
         case (SC_RHOE_VAL); sensor % sVar = SC_RHOE_ID
         case (SC_U_VAL);    sensor % sVar = SC_U_ID
         case (SC_V_VAL);    sensor % sVar = SC_V_ID
         case (SC_W_VAL);    sensor % sVar = SC_W_ID
         case (SC_P_VAL);    sensor % sVar = SC_P_ID
         case (SC_RHOP_VAL); sensor % sVar = SC_RHOP_ID
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
            errorMessage(STD_OUT)
            stop
         end select

      else
         sensor % sVar = SC_RHOP_ID

      end if
!
!     Sensor parameters
!     -----------------
      sensor % s0   = (sensor % high + sensor % low) / 2.0_RP
      sensor % ds   = (sensor % high - sensor % low)
      sensor % ds2  = sensor % ds / 2.0_RP

   end subroutine Set_SCsensor
!
!///////////////////////////////////////////////////////////////////////////////
!
   subroutine Construct_TEsensor(sensor, controlVariables, sem, complementaryGeometry, &
                                 ComputeTimeDerivative, ComputeTimeDerivativeIsolated)
!
!     -------
!     Modules
!     -------
      use FTValueDictionaryClass
      use TransfiniteMapClass, only: TransfiniteHexMap
      use NodalStorageClass,   only: NodalStorage
!
!     ---------
!     Interface
!     ---------
      implicit none
      type(SCsensor_t),                   intent(inout) :: sensor
      type(FTValueDictionary),            intent(in)    :: controlVariables
      type(DGSem),                        intent(in)    :: sem
      logical,                            intent(in)    :: complementaryGeometry
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
      type(TransfiniteHexMap), pointer :: hexMap    => null()
      type(TransfiniteHexMap), pointer :: hex8Map   => null()
      type(TransfiniteHexMap), pointer :: genHexMap => null()

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
!     Generate the complementary geometry if needed
!     ---------------------------------------------
      if (complementaryGeometry) then
      associate(mesh => sensor % TEestim % coarseSem % mesh)

         allocate(hex8Map)
         call hex8Map % constructWithCorners(mesh % elements(1) % SurfInfo % corners)
         allocate(genHexMap)

         do eID = 1, mesh % no_of_elements

            associate(e => mesh % elements(eID))

            if (e % SurfInfo % IsHex8) then
               call hex8Map % setCorners(e % SurfInfo % corners)
               hexMap => hex8Map
            else
               call genHexMap % destruct()
               call genHexMap % constructWithFaces(e % SurfInfo % facePatches)
               hexMap => genHexMap
            end if

            call e % geom % updateComplementaryGrid(NodalStorage(e % Nxyz(IX)), &
                                                    NodalStorage(e % Nxyz(IY)), &
                                                    NodalStorage(e % Nxyz(IZ)), &
                                                    hexMap)

            end associate

         end do

         deallocate(hex8Map)
         deallocate(genHexMap)
         nullify(hexMap)

      end associate
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
         case (SC_ZERO_ID);    write(STD_OUT,"(A)") SC_ZERO_VAL
         case (SC_ONE_ID);     write(STD_OUT,"(A)") SC_ONE_VAL
         case (SC_GRADRHO_ID); write(STD_OUT,"(A)") SC_GRADRHO_VAL
         case (SC_MODAL_ID);   write(STD_OUT,"(A)") SC_MODAL_VAL
         case (SC_TE_ID);      write(STD_OUT,"(A)") SC_TE_VAL
      end select

      if (sensor % sens_type == SC_TE_ID) then
         write(STD_OUT,"(30X,A,A30)", advance="no") "->", "Time derivative: "
         select case (sensor % TEestim % derivType)
            case (ISOLATED_TE);     write(STD_OUT,"(A)") SC_ISOLATED_KEY
            case (NON_ISOLATED_TE); write(STD_OUT,"(A)") SC_NON_ISOLATED_KEY
         end select

         write(STD_OUT,"(30X,A,A30,I0,A,I0,A)") "->", "Delta N: ", sensor % TEestim % deltaN, &
                                                " (min. order is ", sensor % TEestim % Nmin, ")"

      else if (sensor % sens_type /= SC_GRADRHO_ID) then
         write(STD_OUT,"(30X,A,A30)", advance="no") "->", "Sensed variable: "
         select case (sensor % sVar)
            case (SC_RHO_ID);  write(STD_OUT,"(A)") SC_RHO_VAL
            case (SC_RHOU_ID); write(STD_OUT,"(A)") SC_RHOU_VAL
            case (SC_RHOV_ID); write(STD_OUT,"(A)") SC_RHOV_VAL
            case (SC_RHOW_ID); write(STD_OUT,"(A)") SC_RHOW_VAL
            case (SC_RHOE_ID); write(STD_OUT,"(A)") SC_RHOE_VAL
            case (SC_U_ID);    write(STD_OUT,"(A)") SC_U_VAL
            case (SC_V_ID);    write(STD_OUT,"(A)") SC_V_VAL
            case (SC_W_ID);    write(STD_OUT,"(A)") SC_W_VAL
            case (SC_P_ID);    write(STD_OUT,"(A)") SC_P_VAL
            case (SC_RHOP_ID); write(STD_OUT,"(A)") SC_RHOP_VAL
         end select

      end if

      write(STD_OUT,"(30X,A,A30,F5.1)") "->", "Minimum value: ", sensor % low
      write(STD_OUT,"(30X,A,A30,F5.1)") "->", "Maximum value: ", sensor % high


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


      if (allocated(sensor % TEestim))  deallocate(sensor % TEestim)
      if (associated(sensor % Compute)) nullify(sensor % Compute)
      if (associated(sensor % Rescale)) nullify(sensor % Rescale)

   end subroutine Destruct_SCsensor
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


!$omp parallel do schedule(runtime)
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


!$omp parallel do schedule(runtime)
      do eID = 1, sem % mesh % no_of_elements
         sem % mesh % elements(eID) % storage % sensor = 1.0_RP
      end do
!$omp end parallel do

   end subroutine Sensor_one
!
!///////////////////////////////////////////////////////////////////////////////
!
   subroutine Sensor_rho(sensor, sem, t)
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
      real(RP) :: grad_rho2
      real(RP) :: val

!$omp parallel do schedule(runtime)
      do eID = 1, sem % mesh % no_of_elements
      associate(e => sem % mesh % elements(eID))

         val = 0.0_RP
         do k = 0, e % Nxyz(3) ; do j = 0, e % Nxyz(2) ; do i = 0, e % Nxyz(1)
!
!           Compute the square of the norm of the density gradient
!           ------------------------------------------------------
            if ( grad_vars == GRADVARS_STATE ) then
               grad_rho2 = POW2(e % storage % U_x(1,i,j,k)) &
                         + POW2(e % storage % U_y(1,i,j,k)) &
                         + POW2(e % storage % U_z(1,i,j,k))
            else if (grad_vars == GRADVARS_ENTROPY) then
               grad_rho2 = POW2(sum(e % storage % Q(:,i,j,k) * e % storage % U_x(:,i,j,k))) &
                         + POW2(sum(e % storage % Q(:,i,j,k) * e % storage % U_y(:,i,j,k))) &
                         + POW2(sum(e % storage % Q(:,i,j,k) * e % storage % U_z(:,i,j,k)))
            end if
!
!           Integral of the squared gradient
!           --------------------------------
            val = val + NodalStorage(e % Nxyz(1)) % w(i) &
                      * NodalStorage(e % Nxyz(2)) % w(j) &
                      * NodalStorage(e % Nxyz(3)) % w(k) &
                      * e % geom % jacobian(i,j,k)       &
                      * grad_rho2

         end do                ; end do                ; end do

         e % storage % sensor = sqrt(val)

      end associate
      end do
!$omp end parallel do

   end subroutine Sensor_rho
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
      integer               :: Nx, Ny, Nz
      real(RP), allocatable :: sVar(:,:,:)
      real(RP), allocatable :: sVarMod(:,:,:)


      ! Only allocate once
      Nx = maxval(sem % mesh % Nx)
      Ny = maxval(sem % mesh % Ny)
      Nz = maxval(sem % mesh % Nz)
      allocate(sVar(Nx,Ny,Nz))
      allocate(sVarMod(Nx,Ny,Nz))

!$omp parallel do schedule(runtime) private(eID, sVar, sVarMod)
      do eID = 1, sem % mesh % no_of_elements
      associate(e => sem % mesh % elements(eID))

         associate(Nx => e % Nxyz(1), Ny => e % Nxyz(2), Nz => e % Nxyz(3))
!
!        Compute the sensed variable
!        ---------------------------
         do k = 0, Nz ; do j = 0, Ny ; do i = 0, Nx
            sVar(i,j,k) = GetSensedVariable(sensor % sVar, e % storage % Q(:,i,j,k))
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
         if (AlmostEqual(sVarMod(Nx,Ny,Nz), 0.0_RP)) then
            e % storage % sensor = -999.0_RP  ! This is likely to be big enough ;)
         else
            e % storage % sensor = log10( sVarMod(Nx,Ny,Nz)**2 / sum(sVarMod**2) )
         end if

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
      type(Element), pointer               :: e  => null()
      type(Element), pointer               :: ce => null()
      procedure(UserDefinedSourceTermNS_f) :: UserDefinedSourceTermNS


      associate(csem => sensor % TEestim % coarseSem)
!
!     Time derivative
!     ---------------
!$omp parallel do schedule(runtime) private(eID, e, ce)
      do eID = 1, sem % mesh % no_of_elements

         ! TODO: ifort does not like associates with OpenMP
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
      do eID = 1, csem % mesh % no_of_elements

         ! TODO: ifort does not like associates with OpenMP
         e  => sem % mesh % elements(eID)
         ce => csem % mesh % elements(eID)

         ! Use G_NS as temporary storage for Qdot
         call Interp3DArrays(NCONS, e % Nxyz, e % storage % QDot, ce % Nxyz, ce % storage % G_NS)

         mTE = 0.0
         do k = 0, ce % Nxyz(IZ) ; do j = 0, ce % Nxyz(IY) ; do i = 0, ce % Nxyz(IX)

            call UserDefinedSourceTermNS(ce % geom % x, ce % storage % Q, t, S, &
                                         thermodynamics, dimensionless, refValues)

            wx = NodalStorage(ce % Nxyz(IX)) % w(i)
            wy = NodalStorage(ce % Nxyz(IY)) % w(i)
            wz = NodalStorage(ce % Nxyz(IZ)) % w(i)
            Jac = ce % geom % jacobian(i,j,k)

            do iEq = 1, NCONS
               mTE = max(mTE, Jac*wx*wy*wz * abs(ce % storage % QDot(iEq,i,j,k) + S(iEq) - &
                                                 ce % storage % G_NS(iEq,i,j,k)))
            end do

         end do                  ; end do                  ; end do

         if (AlmostEqual(mTE, 0.0_RP)) then
            e % storage % sensor = -999.0_RP  ! This is likely to be big enough ;)
         else
            e % storage % sensor = log10(mTE)
         end if

      end do
!$omp end parallel do

      if (associated(e))  nullify(e)
      if (associated(ce)) nullify(ce)

      end associate

   end subroutine Sensor_truncation
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
   pure function GetSensedVariable(varType, Q) result(s)
!
!     -------
!     Modules
!     -------
      use PhysicsStorage_NS,     only: IRHO, IRHOU, IRHOV, IRHOW, IRHOE
      use VariableConversion_NS, only: Pressure
      implicit none
!
!     ---------
!     Interface
!     ---------
      integer,  intent(in) :: varType
      real(RP), intent(in) :: Q(:)
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

      end select

   end function GetSensedVariable

end module SCsensorClass
