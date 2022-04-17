!
!//////////////////////////////////////////////////////
!
!   @File:    ShockCapturing.f90
!   @Author:  Andrés Mateo (andres.mgabin@upm.es)
!   @Created: Thu Jun 17 2021
!   @Last revision date: Thu Jun 17 2021
!   @Last revision author: Andrés Mateo (andres.mgabin@upm.es)
!
!//////////////////////////////////////////////////////
!
#include "Includes.h"
module ShockCapturing

   use SMConstants,                only: RP, STD_OUT, LINE_LENGTH, NDIM, IX, IY, IZ
   use PhysicsStorage,             only: NCONS, NGRAD
   use FluidData,                  only: dimensionless
   use Utilities,                  only: toLower
   use DGSEMClass,                 only: ComputeTimeDerivative_f, DGSem
   use HexMeshClass,               only: HexMesh
   use ElementClass,               only: Element
#if !defined(SPALARTALMARAS)
   use LESModels,                  only: Smagorinsky_t
#endif
   use SpectralVanishingViscosity, only: SVV, InitializeSVV
   use DGIntegrals,                only: ScalarWeakIntegrals, ScalarStrongIntegrals

   use ShockCapturingKeywords
   use SCsensorClass, only: SCsensor_t, Set_SCsensor, Destruct_SCsensor

   implicit none

   public :: Initialize_ShockCapturing
   public :: ShockCapturingDriver

   type SC_generalMethod_t

         integer, private :: region  !< Sensor region where it acts (1 or 2)

      contains

         procedure :: Initialize => Method_initialize
         procedure :: Viscosity  => Method_viscosity
         procedure :: Advection  => Method_advection
         procedure :: Describe   => Method_describe

   end type SC_generalMethod_t

   type SCdriver_t

      logical  :: isActive = .false.  !< On/Off flag

      type(SCsensor_t),          private              :: sensor   !< Sensor to find discontinuities
      class(SC_generalMethod_t), private, allocatable :: method1  !< First SC method
      class(SC_generalMethod_t), private, allocatable :: method2  !< Second SC method

      contains

         procedure :: Detect           => SC_detect
         procedure :: ComputeViscosity => SC_viscosity
         procedure :: ComputeAdvection => SC_advection
         procedure :: Describe         => SC_describe

         final :: SC_destruct

   end type SCdriver_t

   type, extends(SC_generalMethod_t) :: ArtificialViscosity_t

      integer,  private :: updateMethod     !< Method to compute the viscosity coefficient
      real(RP), private :: mu1              !< First viscosity parameter (low)
      real(RP), private :: alpha1           !< Second viscosity parameter (low)
      real(RP), private :: mu2              !< First viscosity parameter (high)
      real(RP), private :: alpha2           !< Second viscosity parameter (high)
      real(RP), private :: mu2alpha         !< Ratio alpha/mu
      logical,  private :: alphaIsPropToMu  !< .true. if alpha/mu is defined
#if !defined (SPALARTALMARAS)
      type(Smagorinsky_t), private :: Smagorinsky  !< For automatic viscosity
#endif

      contains

         procedure :: Initialize => AV_initialize

   end type ArtificialViscosity_t

   type, extends(SC_generalMethod_t) :: ArtificialAdvection_t
   end type ArtificialAdvection_t

   type, extends(ArtificialViscosity_t) :: SC_NoSVV_t

      integer,                                 private :: fluxType
      procedure(Viscous_Int), nopass, pointer, private :: ViscousFlux => null()

      contains

         procedure :: Initialize => NoSVV_initialize
         procedure :: Viscosity  => NoSVV_viscosity
         procedure :: Describe   => NoSVV_describe

         final :: NoSVV_destruct

   end type SC_NoSVV_t

   type, extends(ArtificialViscosity_t) :: SC_SVV_t

      real(RP), private :: sqrt_mu1
      real(RP), private :: sqrt_alpha1
      real(RP), private :: sqrt_mu2
      real(RP), private :: sqrt_alpha2
      real(RP), private :: sqrt_mu2alpha

      contains

         procedure :: Initialize => SVV_initialize
         procedure :: Viscosity  => SVV_viscosity
         procedure :: Describe   => SVV_describe

   end type SC_SVV_t

   type, extends(ArtificialAdvection_t) :: SC_SSFV_t

      real(RP), private :: c1  ! Lower limit of the blending parameter ("More FV")
      real(RP), private :: c2  ! Higher limit of the blending parameter ("More FS")

      contains

         procedure :: Initialize => SSFV_initialize
         procedure :: Advection  => SSFV_hyperbolic
         procedure :: Describe   => SSFV_describe

   end type SC_SSFV_t
!
!  Interfaces
!  ----------
   abstract interface
      pure subroutine Viscous_Int(nEqn, nGradEqn, Q, Q_x, Q_y, Q_z, mu, beta, kappa, F)
         import RP, NDIM
         integer,       intent(in)  :: nEqn
         integer,       intent(in)  :: nGradEqn
         real(kind=RP), intent(in)  :: Q   (1:nEqn    )
         real(kind=RP), intent(in)  :: Q_x (1:nGradEqn)
         real(kind=RP), intent(in)  :: Q_y (1:nGradEqn)
         real(kind=RP), intent(in)  :: Q_z (1:nGradEqn)
         real(kind=RP), intent(in)  :: mu
         real(kind=RP), intent(in)  :: beta
         real(kind=RP), intent(in)  :: kappa
         real(kind=RP), intent(out) :: F(1:nEqn, 1:NDIM)
      end subroutine Viscous_Int
   end interface

   type(SCdriver_t), allocatable :: ShockCapturingDriver
!
!  ========
   contains
!  ========
!
!
!///////////////////////////////////////////////////////////////////////////////
!
!     Initializer
!
!///////////////////////////////////////////////////////////////////////////////
!
   subroutine Initialize_ShockCapturing(self, controlVariables, sem, &
                                        TimeDerivative, TimeDerivativeIsolated)
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
      type(SCdriver_t), allocatable, intent(inout) :: self
      class(FTValueDictionary),      intent(in)    :: controlVariables
      class(DGSem),                  intent(inout) :: sem
      procedure(ComputeTimeDerivative_f)           :: TimeDerivative
      procedure(ComputeTimeDerivative_f)           :: TimeDerivativeIsolated
!
!     ---------------
!     Local variables
!     ---------------
      character(len=:), allocatable :: method
      logical                       :: updateCompGeometry

!
!     Check if shock-capturing is requested
!     -------------------------------------
      allocate(self)

      if (controlVariables % containsKey(SC_KEY)) then
         self % isActive = controlVariables % logicalValueForKey(SC_KEY)
      else
         self % isActive = .false.
      end if

      if (.not. self % isActive) then
         return
      end if
!
!     Shock-capturing methods
!     -----------------------
      if (controlVariables % containsKey(SC_METHOD1_KEY)) then
         method = controlVariables % stringValueForKey(SC_METHOD1_KEY, LINE_LENGTH)
      else
         method = SC_NO_VAL
      end if
      call toLower(method)

      select case (trim(method))
      case (SC_NOSVV_VAL)
         allocate(SC_NoSVV_t :: self % method1)
         updateCompGeometry = .false.

      case (SC_SVV_VAL)
         allocate(SC_SVV_t :: self % method1)
         updateCompGeometry = .false.

      case (SC_SSFV_VAL)
         allocate(SC_SSFV_t :: self % method1)
         updateCompGeometry = .true.

      case (SC_NO_VAL)
         updateCompGeometry = .false.

      case default
         write(STD_OUT,*) 'ERROR. Unavailable first shock-capturing method. Options are:'
         write(STD_OUT,*) '   * ', SC_NO_VAL
         write(STD_OUT,*) '   * ', SC_NOSVV_VAL
         write(STD_OUT,*) '   * ', SC_SVV_VAL
         write(STD_OUT,*) '   * ', SC_SSFV_VAL
         stop

      end select

      if (controlVariables % containsKey(SC_METHOD2_KEY)) then
         method = controlVariables % stringValueForKey(SC_METHOD2_KEY, LINE_LENGTH)
      else
         method = SC_NO_VAL
      end if
      call toLower(method)

      select case (trim(method))
      case (SC_NOSVV_VAL)
         allocate(SC_NoSVV_t :: self % method2)
         updateCompGeometry = updateCompGeometry .or. .false.

      case (SC_SVV_VAL)
         allocate(SC_SVV_t :: self % method2)
         updateCompGeometry = updateCompGeometry .or. .false.

      case (SC_SSFV_VAL)
         allocate(SC_SSFV_t :: self % method2)
         updateCompGeometry = .true.

      case (SC_NO_VAL)
         updateCompGeometry = updateCompGeometry .or. .false.

      case default
         write(STD_OUT,*) 'ERROR. Unavailable second shock-capturing method. Options are:'
         write(STD_OUT,*) '   * ', SC_NO_VAL
         write(STD_OUT,*) '   * ', SC_NOSVV_VAL
         write(STD_OUT,*) '   * ', SC_SVV_VAL
         write(STD_OUT,*) '   * ', SC_SSFV_VAL
         stop

      end select
!
!     Initialize viscous and hyperbolic terms
!     ---------------------------------------
      if (allocated(self % method1)) call self % method1 % Initialize(controlVariables, sem % mesh, 1)
      if (allocated(self % method2)) call self % method2 % Initialize(controlVariables, sem % mesh, 2)
!
!     Construct sensor
!     ----------------
      call Set_SCsensor(self % sensor, controlVariables, sem, updateCompGeometry, &
                        TimeDerivative, TimeDerivativeIsolated)

   end subroutine Initialize_ShockCapturing
!
!///////////////////////////////////////////////////////////////////////////////
!
!     Base class
!
!///////////////////////////////////////////////////////////////////////////////
!
   subroutine SC_detect(self, sem, t)
!
!     ---------
!     Interface
!     ---------
      implicit none
      class(SCdriver_t), intent(inout) :: self
      type(DGSem),       intent(inout) :: sem
      real(RP),          intent(in)    :: t


      call self % sensor % Compute(sem, t)

   end subroutine SC_detect
!
!///////////////////////////////////////////////////////////////////////////////
!
   subroutine SC_viscosity(self, mesh, e, SCflux)
!
!     ---------
!     Interface
!     ---------
      implicit none
      class(SCdriver_t), intent(in)    :: self
      type(HexMesh),     intent(inout) :: mesh
      type(Element),     intent(inout) :: e
      real(RP),          intent(out)   :: SCflux(1:NCONS,       &
                                                 0:e % Nxyz(1), &
                                                 0:e % Nxyz(2), &
                                                 0:e % Nxyz(3), &
                                                 1:NDIM)
!
!     ---------------
!     Local variables
!     ---------------
      real(RP) :: switch


      switch = self % sensor % Rescale(e % storage % sensor)

      if (switch >= 1.0_RP) then
         if (allocated(self % method2)) then
            call self % method2 % Viscosity(mesh, e, switch, SCflux)
         end if

      else if (switch > 0.0_RP) then
         if (allocated(self % method1)) then
            call self % method1 % Viscosity(mesh, e, switch, SCflux)
         end if

      else
         SCflux = 0.0_RP
         e % storage % artificialDiss = 0.0_RP

      end if

   end subroutine SC_viscosity
!
!///////////////////////////////////////////////////////////////////////////////
!
   function SC_advection(self, e, Fv, Qdot, isStrong) result(computed)
!
!     ---------
!     Interface
!     ---------
      implicit none
      class(SCdriver_t), intent(in)  :: self
      type(Element),     intent(in)  :: e
      logical, optional, intent(in)  :: isStrong
      real(RP),          intent(in)  :: Fv(1:NCONS,       &
                                           0:e % Nxyz(1), &
                                           0:e % Nxyz(2), &
                                           0:e % Nxyz(3), &
                                           1:NDIM         )
      real(RP),          intent(out) :: Qdot(1:NCONS,       &
                                             0:e % Nxyz(1), &
                                             0:e % Nxyz(2), &
                                             0:e % Nxyz(3)  )
      logical :: computed
!
!     ---------------
!     Local variables
!     ---------------
      real(RP) :: switch
      logical  :: strongIntegral


      if (present(isStrong)) then
         strongIntegral = isStrong
      else
         strongIntegral = .false.
      end if

      switch = self % sensor % Rescale(e % storage % sensor)

      if (switch >= 1.0_RP) then
         if (allocated(self % method2)) then
            computed = self % method2 % Advection(e, switch, Fv, strongIntegral, Qdot)
         end if

      else if (switch > 0.0_RP) then
         if (allocated(self % method1)) then
            computed = self % method1 % Advection(e, switch, Fv, strongIntegral, Qdot)
         end if

      else
         computed = .false.

      end if

   end function SC_advection
!
!///////////////////////////////////////////////////////////////////////////////
!
   subroutine SC_describe(self)
!
!     -------
!     Modules
!     -------
      use MPI_Process_Info, only: MPI_Process
      use Headers,          only: Subsection_Header
!
!     ---------
!     Interface
!     ---------
      implicit none
      class(SCdriver_t), intent(in) :: self


      if (.not. MPI_Process % isRoot .or. .not. self % isActive) return

      write(STD_OUT, "(/)")
      call Subsection_Header("Shock-Capturing")

      call self % sensor % Describe()

      if (allocated(self % method1)) then
         write(STD_OUT,*) ""
         write(STD_OUT,"(30X,A,A30)") "=>", "First method"
         call self % method1 % Describe()
      end if

      if (allocated(self % method2)) then
         write(STD_OUT,*) ""
         write(STD_OUT,"(30X,A,A30)") "=>", "Second method"
         call self % method2 % Describe()
      end if

   end subroutine SC_describe
!
!///////////////////////////////////////////////////////////////////////////////
!
   subroutine SC_destruct(self)
!
!     ---------
!     Interface
!     ---------
      implicit none
      type(SCdriver_t), intent(inout) :: self


      self%isActive = .false.

      call Destruct_SCsensor(self % sensor)

      if (allocated(self % method1)) deallocate(self % method1)
      if (allocated(self % method2)) deallocate(self % method2)

   end subroutine SC_destruct
!
!///////////////////////////////////////////////////////////////////////////////
!
!     Shock-capturing methods base class
!
!///////////////////////////////////////////////////////////////////////////////
!
   subroutine Method_initialize(self, controlVariables, mesh, region)
!
!     -------
!     Modules
!     -------
      use FTValueDictionaryClass
!
!     ---------
!     Interface
!     ---------
      class(SC_generalMethod_t), intent(inout) :: self
      type(FTValueDictionary),   intent(in)    :: controlVariables
      type(HexMesh),             intent(inout) :: mesh
      integer,                   intent(in)    :: region

   end subroutine Method_initialize
!
!///////////////////////////////////////////////////////////////////////////////
!
   subroutine Method_viscosity(self, mesh, e, switch, SCflux)
!
!     ---------
!     Interface
!     ---------
      implicit none
      class(SC_generalMethod_t), intent(in)    :: self
      type(HexMesh),             intent(inout) :: mesh
      type(Element),             intent(inout) :: e
      real(RP),                  intent(in)    :: switch
      real(RP),                  intent(out)   :: SCflux(1:NCONS,       &
                                                         0:e % Nxyz(1), &
                                                         0:e % Nxyz(2), &
                                                         0:e % Nxyz(3), &
                                                         1:NDIM)


      SCflux = 0.0_RP

   end subroutine Method_viscosity
!
!///////////////////////////////////////////////////////////////////////////////
!
   function Method_advection(self, e, switch, Fv, isStrong, div)
!
!     ---------
!     Interface
!     ---------
      implicit none
      class(SC_generalMethod_t), intent(in)    :: self
      type(Element),             intent(in)    :: e
      real(RP),                  intent(in)    :: switch
      logical,                   intent(in)    :: isStrong
      real(RP),                  intent(in)    :: Fv(1:NCONS,       &
                                                     0:e % Nxyz(1), &
                                                     0:e % Nxyz(2), &
                                                     0:e % Nxyz(3), &
                                                     1:NDIM         )
      real(RP),                  intent(inout) :: div(1:NCONS,       &
                                                      0:e % Nxyz(1), &
                                                      0:e % Nxyz(2), &
                                                      0:e % Nxyz(3)  )
      logical :: Method_advection


      Method_advection = .false.

   end function Method_advection
!
!///////////////////////////////////////////////////////////////////////////////
!
   subroutine Method_describe(self)
!
!     ---------
!     Interface
!     ---------
      implicit none
      class(SC_generalMethod_t), intent(in) :: self

   end subroutine Method_describe
!
!///////////////////////////////////////////////////////////////////////////////
!
   subroutine AV_initialize(self, controlVariables, mesh, region)
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
      class(ArtificialViscosity_t), intent(inout) :: self
      type(FTValueDictionary),      intent(in)    :: controlVariables
      type(HexMesh),                intent(inout) :: mesh
      integer,                      intent(in)    :: region
!
!     ---------------
!     Local variables
!     ---------------
      character(len=:), allocatable :: update

!
!     Viscosity values (mu and alpha)
!     -------------------------------
      if (region == 1) then

         if (controlVariables % containsKey(SC_MU1_KEY)) then
            self % mu1 = controlVariables % doublePrecisionValueForKey(SC_MU1_KEY)
         else
            write(STD_OUT,*) "ERROR. A value for the artificial 'mu 1' must be given."
            stop
         end if

         if (controlVariables % containsKey(SC_MU2_KEY)) then
            self % mu2 = controlVariables % doublePrecisionValueForKey(SC_MU2_KEY)
         else
            self % mu2 = self % mu1
         end if

      else

         if (controlVariables % containsKey(SC_MU2_KEY)) then
            self % mu2 = controlVariables % doublePrecisionValueForKey(SC_MU2_KEY)
         else
            write(STD_OUT,*) "ERROR. A value for the artificial 'mu 2' must be given."
            stop
         end if
         self % mu1 = self % mu2

      end if

      if (controlVariables % containsKey(SC_ALPHA_MU_KEY)) then

         self % alphaIsPropToMu = .true.
         self % mu2alpha        = controlVariables % doublePrecisionValueForKey(SC_ALPHA_MU_KEY)
         self % alpha1          = self % mu2alpha * self % mu1
         self % alpha2          = self % mu2alpha * self % mu2

      else

         self % alphaIsPropToMu = .false.

         if (region == 1) then

            if (controlVariables % containsKey(SC_ALPHA1_KEY)) then
               self % alpha1 = controlVariables % doublePrecisionValueForKey(SC_ALPHA1_KEY)
            else
               self % alpha1 = 0.0_RP
            end if

            if (controlVariables % containsKey(SC_ALPHA2_KEY)) then
               self % alpha2 = controlVariables % doublePrecisionValueForKey(SC_ALPHA2_KEY)
            else
               self % alpha2 = self % alpha1
            end if

         else

            if (controlVariables % containsKey(SC_ALPHA2_KEY)) then
               self % alpha2 = controlVariables % doublePrecisionValueForKey(SC_ALPHA2_KEY)
            else
               self % alpha2 = 0.0_RP
            end if
            self % alpha1 = self % alpha2

         end if

      end if
!
!     Viscosity update method
!     -----------------------
      if (region == 1) then

         if (controlVariables % containsKey(SC_UPDATE_KEY)) then
            update = controlVariables % StringValueForKey(SC_UPDATE_KEY, LINE_LENGTH)
            call toLower(update)

            select case (trim(update))
            case (SC_CONST_VAL)
               self % updateMethod = SC_CONST_ID

            case (SC_SENSOR_VAL)
               self % updateMethod = SC_SENSOR_ID

#if !defined (SPALARTALMARAS)
            case (SC_SMAG_VAL)

               self % updateMethod = SC_SMAG_ID
               if (.not. self%alphaIsPropToMu) then
                  write(STD_OUT,*) 'ERROR. Alpha must be proportional to mu when using shock-capturing with LES.'
                  stop
               end if

               ! TODO: Use the default constructor
               self % Smagorinsky % active = .true.
               self % Smagorinsky % requiresWallDistances = .false.
               self % Smagorinsky % WallModel = 0  ! No wall model
               self % Smagorinsky % CS = self % mu1
#endif

            case default
               write(STD_OUT,*) 'ERROR. Unavailable shock-capturing update strategy. Options are:'
               write(STD_OUT,*) '   * ', SC_CONST_VAL
               write(STD_OUT,*) '   * ', SC_SENSOR_VAL
#if !defined (SPALARTALMARAS)
               write(STD_OUT,*) '   * ', SC_SMAG_VAL
#endif
               stop

            end select

         else
            self % updateMethod = SC_CONST_ID

         end if

      else

         self % updateMethod = SC_CONST_ID

      end if

   end subroutine AV_initialize
!
!///////////////////////////////////////////////////////////////////////////////
!
!     Simple artificial viscosity
!
!///////////////////////////////////////////////////////////////////////////////
!
   subroutine NoSVV_initialize(self, controlVariables, mesh, region)
!
!     -------
!     Modules
!     -------
      use FTValueDictionaryClass
      use PhysicsStorage, only: grad_vars, GRADVARS_STATE, &
                                GRADVARS_ENTROPY, GRADVARS_ENERGY
      use Physics,        only: ViscousFlux_STATE, ViscousFlux_ENTROPY, &
                                ViscousFlux_ENERGY, GuermondPopovFlux_ENTROPY
!
!     ---------
!     Interface
!     ---------
      implicit none
      class(SC_NoSVV_t),       intent(inout) :: self
      type(FTValueDictionary), intent(in)    :: controlVariables
      type(HexMesh),           intent(inout) :: mesh
      integer,                 intent(in)    :: region
!
!     ---------------
!     Local variables
!     ---------------
      character(len=:), allocatable :: flux

!
!     Parent initializer
!     ------------------
      call self % ArtificialViscosity_t % Initialize(controlVariables, mesh, region)
      self % region = region
!
!     Set the flux type
!     -----------------
      if (region == 1 .and. controlVariables % containsKey(SC_VISC_FLUX1_KEY)) then
         flux = controlVariables % stringValueForKey(SC_VISC_FLUX1_KEY, LINE_LENGTH)
      else if (region == 2 .and. controlVariables % containsKey(SC_VISC_FLUX2_KEY)) then
         flux = controlVariables % stringValueForKey(SC_VISC_FLUX2_KEY, LINE_LENGTH)
      else
         flux = SC_PHYS_VAL
      end if

      call toLower(flux)

      select case (trim(flux))
      case (SC_PHYS_VAL); self % fluxType = SC_PHYS_ID
      case (SC_GP_VAL);   self % fluxType = SC_GP_ID
      case default
         write(STD_OUT,'(A,I1,A)') 'ERROR. Artificial viscosity ', region, &
                                   ' not recognized. Options are:'
         write(STD_OUT,*) '   * ', SC_PHYS_VAL
         write(STD_OUT,*) '   * ', SC_GP_VAL
         stop
      end select

      select case (self % fluxType)
      case (SC_PHYS_ID)
         select case (grad_vars)
         case (GRADVARS_STATE);   self % ViscousFlux => ViscousFlux_STATE
         case (GRADVARS_ENTROPY); self % ViscousFlux => ViscousFlux_ENTROPY
         case (GRADVARS_ENERGY);  self % ViscousFlux => ViscousFlux_ENERGY
         end select

      case (SC_GP_ID)
         select case (grad_vars)
         case (GRADVARS_ENTROPY); self % ViscousFlux => GuermondPopovFlux_ENTROPY
         case default
            write(STD_OUT,*) "ERROR. Guermond-Popov (2014) artificial ",  &
                              "viscosity is only configured for Entropy ", &
                              "gradient variables"
            stop
         end select

      end select

   end subroutine NoSVV_initialize
!
!///////////////////////////////////////////////////////////////////////////////
!
   subroutine NoSVV_viscosity(self, mesh, e, switch, SCflux)
!
!     --------------------------------------------------------------------------
!     TODO: Introduce alpha viscosity, which probably means reimplementing here
!           all the viscous fluxes of `Physics_NS`...
!     --------------------------------------------------------------------------
!
!     ---------
!     Interface
!     ---------
      implicit none
      class(SC_NoSVV_t), intent(in)    :: self
      type(HexMesh),     intent(inout) :: mesh
      type(Element),     intent(inout) :: e
      real(RP),          intent(in)    :: switch
      real(RP),          intent(out)   :: SCflux(1:NCONS,       &
                                                 0:e % Nxyz(1), &
                                                 0:e % Nxyz(2), &
                                                 0:e % Nxyz(3), &
                                                 1:NDIM)
!
!     ---------------
!     Local variables
!     ---------------
      integer  :: i
      integer  :: j
      integer  :: k
      integer  :: fIDs(6)
      real(RP) :: delta
      real(RP) :: kappa
      real(RP) :: mu(0:e % Nxyz(1), 0:e % Nxyz(2), 0:e % Nxyz(3))
      real(RP) :: covariantFlux(1:NCONS, 1:NDIM)


      if (switch > 0.0_RP) then
!
!        Compute viscosity
!        -----------------
         select case (self % updateMethod)
         case (SC_CONST_ID)
            mu = merge(self % mu2, self % mu1, switch >= 1.0_RP)

         case (SC_SENSOR_ID)
            mu = self % mu1 * (1.0_RP-switch) + self % mu2 * switch

#if !defined (SPALARTALMARAS)
         case (SC_SMAG_ID)

            delta = (e % geom % Volume / product(e % Nxyz + 1)) ** (1.0_RP / 3.0_RP)
            do k = 0, e % Nxyz(3) ; do j = 0, e % Nxyz(2) ; do i = 0, e % Nxyz(1)
               call self % Smagorinsky % ComputeViscosity(delta, e % geom % dWall(i,j,k), &
                                                          e % storage % Q(:,i,j,k),       &
                                                          e % storage % U_x(:,i,j,k),     &
                                                          e % storage % U_y(:,i,j,k),     &
                                                          e % storage % U_z(:,i,j,k),     &
                                                          mu(i,j,k))
            end do                ; end do                ; end do
#endif

         end select
!
!        Compute the viscous flux
!        ------------------------
         do k = 0, e % Nxyz(3) ; do j = 0, e % Nxyz(2) ; do i = 0, e % Nxyz(1)

            kappa = dimensionless % mu_to_kappa * mu(i,j,k)
            call self % ViscousFlux(NCONS, NGRAD, e % storage % Q(:,i,j,k), &
                                    e % storage % U_x(:,i,j,k),             &
                                    e % storage % U_y(:,i,j,k),             &
                                    e % storage % U_z(:,i,j,k),             &
                                    mu(i,j,k), 0.0_RP, kappa,               &
                                    covariantflux)

            SCflux(:,i,j,k,IX) = covariantFlux(:,IX) * e % geom % jGradXi(IX,i,j,k) &
                               + covariantFlux(:,IY) * e % geom % jGradXi(IY,i,j,k) &
                               + covariantFlux(:,IZ) * e % geom % jGradXi(IZ,i,j,k)


            SCflux(:,i,j,k,IY) = covariantFlux(:,IX) * e % geom % jGradEta(IX,i,j,k) &
                               + covariantFlux(:,IY) * e % geom % jGradEta(IY,i,j,k) &
                               + covariantFlux(:,IZ) * e % geom % jGradEta(IZ,i,j,k)


            SCflux(:,i,j,k,IZ) = covariantFlux(:,IX) * e % geom % jGradZeta(IX,i,j,k) &
                               + covariantFlux(:,IY) * e % geom % jGradZeta(IY,i,j,k) &
                               + covariantFlux(:,IZ) * e % geom % jGradZeta(IZ,i,j,k)

         end do                ; end do                ; end do

      else

         e % storage % artificialDiss = 0.0_RP
         SCflux = 0.0_RP

      end if
!
!     Project to faces
!     ----------------
      fIDs = e % faceIDs
      call e % ProlongAviscFluxToFaces(NCONS, SCflux, mesh % faces(fIDs(1)), &
                                                      mesh % faces(fIDs(2)), &
                                                      mesh % faces(fIDs(3)), &
                                                      mesh % faces(fIDs(4)), &
                                                      mesh % faces(fIDs(5)), &
                                                      mesh % faces(fIDs(6))  )

   end subroutine NoSVV_viscosity
!
!///////////////////////////////////////////////////////////////////////////////
!
   subroutine NoSVV_describe(self)
!
!     -------
!     Modules
!     -------
      use MPI_Process_Info, only: MPI_Process
!
!     ---------
!     Interface
!     ---------
      implicit none
      class(SC_NoSVV_t), intent(in) :: self


      if (.not. MPI_Process % isRoot) return

      write(STD_OUT,"(30X,A,A30)", advance="no") "->", "Viscosity update method: "
      select case (self % updateMethod)
         case (SC_CONST_ID);  write(STD_OUT,"(A)") SC_CONST_VAL
         case (SC_SENSOR_ID); write(STD_OUT,"(A)") SC_SENSOR_VAL
         case (SC_SMAG_ID);   write(STD_OUT,"(A)") SC_SMAG_VAL
      end select

      if (self % updateMethod == SC_SMAG_ID) then
#if !defined (SPALARTALMARAS)
         write(STD_OUT,"(30X,A,A30,F4.2)") "->", "LES intensity (CS): ", self % Smagorinsky % CS
#endif
      else
         write(STD_OUT,"(30X,A,A30,F10.6)") "->","Mu viscosity 1: ", self % mu1
         write(STD_OUT,"(30X,A,A30,F10.6)") "->","Mu viscosity 2: ", self % mu2
      end if

      if (self % alphaIsPropToMu) then
         write(STD_OUT,"(30X,A,A30,F7.3,A)") "->", "Alpha viscosity: ", self % mu2alpha, "x mu"
      else
         write(STD_OUT,"(30X,A,A30,F10.6)") "->","Alpha viscosity 1: ", self % alpha1
         write(STD_OUT,"(30X,A,A30,F10.6)") "->","Alpha viscosity 2: ", self % alpha2
      end if

      write(STD_OUT,"(30X,A,A30)", advance="no") "->", "Dissipation type: "
      select case (self % fluxType)
         case (SC_PHYS_ID); write(STD_OUT,"(A)") SC_PHYS_VAL
         case (SC_GP_ID);   write(STD_OUT,"(A)") SC_GP_VAL
      end select

   end subroutine NoSVV_describe
!
!///////////////////////////////////////////////////////////////////////////////
!
   pure subroutine NoSVV_destruct(self)
!
!     ---------
!     Interface
!     ---------
      implicit none
      type(SC_NoSVV_t), intent(inout) :: self


      if (associated(self % ViscousFlux)) nullify(self % ViscousFlux)

   end subroutine NoSVV_destruct
!
!///////////////////////////////////////////////////////////////////////////////
!
!     SVV filtered viscosity
!
!///////////////////////////////////////////////////////////////////////////////
!
   subroutine SVV_initialize(self, controlVariables, mesh, region)
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
      class(SC_SVV_t),         intent(inout) :: self
      type(FTValueDictionary), intent(in)    :: controlVariables
      type(HexMesh),           intent(inout) :: mesh
      integer,                 intent(in)    :: region


      ! TODO: Implement it also for region 2, but SVV does not seem very useful there...
      if (region == 2) then
         write(STD_OUT,*) "ERROR. SVV viscosity can be used only in the first region of the sensor."
         stop
      end if
!
!     Parent initializer
!     ------------------
      call self % ArtificialViscosity_t % Initialize(controlVariables, mesh, region)
      self % region = region
!
!     Set the square root of the viscosities
!     --------------------------------------
      self % sqrt_mu1 = sqrt(self % mu1)
      self % sqrt_mu2 = sqrt(self % mu2)

      if (self % alphaIsPropToMu) then
         self % sqrt_mu2alpha = sqrt(self % mu2alpha)
      else
         self % sqrt_alpha1 = sqrt(self % alpha1)
         self % sqrt_alpha2 = sqrt(self % alpha2)
      end if
!
!     Start the SVV module
!     --------------------
      call InitializeSVV(SVV, controlVariables, mesh)

   end subroutine SVV_initialize
!
!///////////////////////////////////////////////////////////////////////////////
!
   subroutine SVV_viscosity(self, mesh, e, switch, SCflux)
!
!     ---------
!     Interface
!     ---------
      implicit none
      class(SC_SVV_t), intent(in)    :: self
      type(HexMesh),   intent(inout) :: mesh
      type(Element),   intent(inout) :: e
      real(RP),        intent(in)    :: switch
      real(RP),        intent(out)   :: SCflux(1:NCONS,       &
                                               0:e % Nxyz(1), &
                                               0:e % Nxyz(2), &
                                               0:e % Nxyz(3), &
                                               1:NDIM)
!
!     ---------------
!     Local variables
!     ---------------
      integer  :: i
      integer  :: j
      integer  :: k
      integer  :: fIDs(6)
      real(RP) :: delta
      real(RP) :: salpha
      real(RP) :: sqrt_mu(0:e % Nxyz(1), 0:e % Nxyz(2), 0:e % Nxyz(3))
      real(RP) :: sqrt_alpha(0:e % Nxyz(1), 0:e % Nxyz(2), 0:e % Nxyz(3))


      if (switch > 0.0_RP) then
!
!        Compute viscosities
!        -------------------
         select case (self % updateMethod)
         case (SC_CONST_ID)
            sqrt_mu = merge(self % sqrt_mu2,    self % sqrt_mu1,    switch >= 1.0_RP)
            salpha  = merge(self % sqrt_alpha2, self % sqrt_alpha1, switch >= 1.0_RP)

         case (SC_SENSOR_ID)
            sqrt_mu = self % sqrt_mu1 * (1.0_RP-switch) + self % sqrt_mu2 * switch
            salpha  = self % sqrt_alpha1 * (1.0_RP-switch) + self % sqrt_alpha2 * switch

#if !defined (SPALARTALMARAS)
         case (SC_SMAG_ID)

            delta = (e % geom % Volume / product(e % Nxyz + 1)) ** (1.0_RP / 3.0_RP)
            do k = 0, e % Nxyz(3) ; do j = 0, e % Nxyz(2) ; do i = 0, e % Nxyz(1)
               call self % Smagorinsky % ComputeViscosity(delta, e % geom % dWall(i,j,k), &
                                                          e % storage % Q(:,i,j,k),       &
                                                          e % storage % U_x(:,i,j,k),     &
                                                          e % storage % U_y(:,i,j,k),     &
                                                          e % storage % U_z(:,i,j,k),     &
                                                          sqrt_mu(i,j,k))
               sqrt_mu(i,j,k) = sqrt(sqrt_mu(i,j,k))
            end do                ; end do                ; end do
#endif

         end select

         if (self % alphaIsPropToMu) then
            sqrt_alpha = self % sqrt_mu2alpha * sqrt_mu
         else
            sqrt_alpha = salpha
         end if
!
!        Compute the viscous flux
!        ------------------------
         call SVV % ComputeInnerFluxes(e, sqrt_mu, sqrt_alpha, SCflux)

      else

         e % storage % artificialDiss = 0.0_RP
         SCflux = 0.0_RP

      end if
!
!     Project to faces
!     ----------------
      fIDs = e % faceIDs
      call e % ProlongAviscFluxToFaces(NCONS, SCflux, mesh % faces(fIDs(1)), &
                                                      mesh % faces(fIDs(2)), &
                                                      mesh % faces(fIDs(3)), &
                                                      mesh % faces(fIDs(4)), &
                                                      mesh % faces(fIDs(5)), &
                                                      mesh % faces(fIDs(6))  )

   end subroutine SVV_viscosity
!
!///////////////////////////////////////////////////////////////////////////////
!
   subroutine SVV_describe(self)
!
!     -------
!     Modules
!     -------
      use MPI_Process_Info, only: MPI_Process
!
!     ---------
!     Interface
!     ---------
      implicit none
      class(SC_SVV_t), intent(in) :: self


      if (.not. MPI_Process % isRoot) return

      write(STD_OUT,"(30X,A,A30)", advance="no") "->", "Viscosity update method: "
      select case (self % updateMethod)
         case (SC_CONST_ID);  write(STD_OUT,"(A)") SC_CONST_VAL
         case (SC_SENSOR_ID); write(STD_OUT,"(A)") SC_SENSOR_VAL
         case (SC_SMAG_ID);   write(STD_OUT,"(A)") SC_SMAG_VAL
      end select

      if (self % updateMethod == SC_SMAG_ID) then
#if !defined (SPALARTALMARAS)
         write(STD_OUT,"(30X,A,A30,F4.2)") "->", "LES intensity (CS): ", self % Smagorinsky % CS
#endif
      else
         write(STD_OUT,"(30X,A,A30,F10.6)") "->","Mu viscosity 1: ", self % mu1
         write(STD_OUT,"(30X,A,A30,F10.6)") "->","Mu viscosity 2: ", self % mu2
      end if

      if (self % alphaIsPropToMu) then
         write(STD_OUT,"(30X,A,A30,F7.3,A)") "->", "Alpha viscosity: ", self % mu2alpha, "x mu"
      else
         write(STD_OUT,"(30X,A,A30,F10.6)") "->","Alpha viscosity 1: ", self % alpha1
         write(STD_OUT,"(30X,A,A30,F10.6)") "->","Alpha viscosity 2: ", self % alpha2
      end if

      call SVV % Describe()

   end subroutine SVV_describe
!
!///////////////////////////////////////////////////////////////////////////////
!
!     Entropy stable finite volumes (SSFV)
!
!///////////////////////////////////////////////////////////////////////////////
!
   subroutine SSFV_initialize(self, controlVariables, mesh, region)
!
!     -------
!     Modules
!     -------
      use FTValueDictionaryClass
      use TransfiniteMapClass,   only: TransfiniteHexMap
      use NodalStorageClass,     only: NodalStorage
!
!     ---------
!     Interface
!     ---------
      implicit none
      class(SC_SSFV_t),        intent(inout) :: self
      type(FTValueDictionary), intent(in)    :: controlVariables
      type(HexMesh),           intent(inout) :: mesh
      integer,                 intent(in)    :: region
!
!     ---------------
!     Local variables
!     ---------------
!
      integer                          :: eID
      type(TransfiniteHexMap), pointer :: hexMap    => null()
      type(TransfiniteHexMap), pointer :: hex8Map   => null()
      type(TransfiniteHexMap), pointer :: genHexMap => null()


      self % region = region
!
!     Blending parameters
!     -------------------
      if (region == 1) then

         if (controlVariables % containsKey(SC_BLENDING1_KEY)) then
            self % c1 = controlVariables % doublePrecisionValueForKey(SC_BLENDING1_KEY)
         else
            write(STD_OUT,*) "ERROR. A value for the blending parameter 1 must be given."
            stop
         end if

         if (controlVariables % containsKey(SC_BLENDING2_KEY)) then
            self % c2 = controlVariables % doublePrecisionValueForKey(SC_BLENDING2_KEY)
         else
            self % c2 = self % c1
         end if

      else

         if (controlVariables % containsKey(SC_BLENDING2_KEY)) then
            self % c2 = controlVariables % doublePrecisionValueForKey(SC_BLENDING2_KEY)
         else
            write(STD_OUT,*) "ERROR. A value for the blending parameter 2 must be given."
            stop
         end if
         self % c1 = self % c2

      end if
!
!     Compute the geometry for the complementary grid
!     -----------------------------------------------
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
!
!     Release the memory
!     ------------------
      deallocate(hex8Map)
      deallocate(genHexMap)
      nullify(hexMap)

   end subroutine SSFV_initialize
!
!///////////////////////////////////////////////////////////////////////////////
!
   function SSFV_hyperbolic(self, e, switch, Fv, isStrong, div) result(computed)
!
!     -------
!     Modules
!     -------
      use PhysicsStorage,            only: IRHO
      use NodalStorageClass,         only: NodalStorage
      use HyperbolicDiscretizations, only: HyperbolicDiscretization
      use Physics,                   only: EulerFlux
      use VariableConversion,        only: Pressure, getEntropyVariables
      ! use RiemannSolvers_NS,         only: RiemannSolver
      use SMConstants,               only: PI
!
!     ---------
!     Interface
!     ---------
      implicit none
      class(SC_SSFV_t),    intent(in) :: self
      type(Element),       intent(in) :: e
      real(RP),            intent(in) :: switch
      logical,             intent(in) :: isStrong
      real(RP),            intent(in) :: Fv(1:NCONS,       &
                                            0:e % Nxyz(1), &
                                            0:e % Nxyz(2), &
                                            0:e % Nxyz(3), &
                                            1:NDIM         )
      real(RP),         intent(inout) :: div(1:NCONS,       &
                                             0:e % Nxyz(1), &
                                             0:e % Nxyz(2), &
                                             0:e % Nxyz(3)  )
      logical :: computed
!
!     ---------------
!     Local variables
!     ---------------
      integer  :: i, j, k, r, s
      real(RP) :: Fc    (NCONS, 0:e % Nxyz(1),   0:e % Nxyz(2), 0:e % Nxyz(3), NDIM)
      real(RP) :: Fsharp(NCONS, 0:e % Nxyz(1),   0:e % Nxyz(1), 0:e % Nxyz(2), 0:e % Nxyz(3))
      real(RP) :: Gsharp(NCONS, 0:e % Nxyz(2),   0:e % Nxyz(1), 0:e % Nxyz(2), 0:e % Nxyz(3))
      real(RP) :: Hsharp(NCONS, 0:e % Nxyz(3),   0:e % Nxyz(1), 0:e % Nxyz(2), 0:e % Nxyz(3))
      real(RP) :: FSx   (NCONS, 0:e % Nxyz(1)+1, 0:e % Nxyz(2), 0:e % Nxyz(3))
      real(RP) :: FSy   (NCONS, 0:e % Nxyz(2)+1, 0:e % Nxyz(1), 0:e % Nxyz(3))
      real(RP) :: FSz   (NCONS, 0:e % Nxyz(3)+1, 0:e % Nxyz(1), 0:e % Nxyz(2))
      real(RP) :: FVx   (NCONS, 1:e % Nxyz(1)  , 0:e % Nxyz(2), 0:e % Nxyz(3))
      real(RP) :: FVy   (NCONS, 1:e % Nxyz(2)  , 0:e % Nxyz(1), 0:e % Nxyz(3))
      real(RP) :: FVz   (NCONS, 1:e % Nxyz(3)  , 0:e % Nxyz(1), 0:e % Nxyz(2))
      real(RP) :: w1(NCONS)
      real(RP) :: w2(NCONS)
      real(RP) :: p
      real(RP) :: invRho
      real(RP) :: b
      real(RP) :: c2
      real(RP) :: d


      associate(Nx => e % Nxyz(1), &
                Ny => e % Nxyz(2), &
                Nz => e % Nxyz(3)  )
      associate(spAxi   => NodalStorage(Nx), &
                spAeta  => NodalStorage(Ny), &
                spAzeta => NodalStorage(Nz)  )
!
!     Entropy-conservative averaging in complementary points
!     ------------------------------------------------------
      call HyperbolicDiscretization % ComputeInnerFluxes(e, EulerFlux, Fc)
      call HyperbolicDiscretization % ComputeSplitFormFluxes(e, Fc, Fsharp, Gsharp, Hsharp)

      do k = 0, Nz ; do j = 0, Ny ; do i = 1, Nx
         FSx(:,i,j,k) = 0.0_RP
         do s = i, Nx ; do r = 0, i-1
            FSx(:,i,j,k) = FSx(:,i,j,k) + spAxi % Q(r,s) * Fsharp(:,r,s,j,k)
         end do                  ; end do
         FSx(:,i,j,k) = 2.0_RP * FSx(:,i,j,k)
      end do                ; end do                ; end do

      do k = 0, Nz ; do i = 0, Nx ; do j = 1, Ny
         FSy(:,j,i,k) = 0.0_RP
         do s = j, Ny ; do r = 0, j-1
            FSy(:,j,i,k) = FSy(:,j,i,k) + spAeta % Q(r,s) * Gsharp(:,r,i,s,k)
         end do                  ; end do
         FSy(:,j,i,k) = 2.0_RP * FSy(:,j,i,k)
      end do                ; end do                ; end do

      do j = 0, Ny ; do i = 0, Nx ; do k = 1, Nz
         FSz(:,k,i,j) = 0.0_RP
         do s = k, Nz ; do r = 0, k-1
            FSz(:,k,i,j) = FSz(:,k,i,j) + spAzeta % Q(r,s) * Hsharp(:,r,i,j,s)
         end do                  ; end do
         FSz(:,k,i,j) = 2.0_RP * FSz(:,k,i,j)
      end do                ; end do                ; end do
      ! do k = 0, Nz ; do j = 0, Ny ; do i = 1, Nx
      !    call RiemannSolver(e % storage % Q(:,i-1,j,k), &
      !                       e % storage % Q(:,i,j,k),   &
      !                       e % geom % ncXi(:,i,j,k),   &
      !                       e % geom % t1cXi(:,i,j,k),  &
      !                       e % geom % t2cXi(:,i,j,k),  &
      !                       FSx(:,i,j,k))
      !    FSx(:,i,j,k) = FSx(:,i,j,k) * e % geom % JfcXi(i,j,k)
      ! end do                ; end do                ; end do

      ! do k = 0, Nz ; do i = 0, Nx ; do j = 1, Ny
      !    call RiemannSolver(e % storage % Q(:,i,j-1,k), &
      !                       e % storage % Q(:,i,j,k),   &
      !                       e % geom % ncEta(:,i,j,k),  &
      !                       e % geom % t1cEta(:,i,j,k), &
      !                       e % geom % t2cEta(:,i,j,k), &
      !                       FSy(:,j,i,k))
      !    FSy(:,j,i,k) = FSy(:,j,i,k) * e % geom % JfcEta(i,j,k)
      ! end do                ; end do                ; end do

      ! do j = 0, Ny ; do i = 0, Nx ; do k = 1, Nz
      !    call RiemannSolver(e % storage % Q(:,i,j,k-1),  &
      !                       e % storage % Q(:,i,j,k),    &
      !                       e % geom % ncZeta(:,i,j,k),  &
      !                       e % geom % t1cZeta(:,i,j,k), &
      !                       e % geom % t2cZeta(:,i,j,k), &
      !                       FSz(:,k,i,j))
      !    FSz(:,k,i,j) = FSz(:,k,i,j) * e % geom % JfcZeta(i,j,k)
      ! end do                ; end do                ; end do
!
!     Boundaries
!     ----------
      do k = 0, Nz ; do j = 0, Ny
         FSx(:,0,j,k) = Fc(:,0,j,k,IX)
         FSx(:,Nx+1,j,k) = Fc(:,Nx,j,k,IX)
      end do       ; end do

      do k = 0, Nz ; do i = 0, Nx
         FSy(:,0,i,k) = Fc(:,i,0,k,IY)
         FSy(:,Ny+1,i,k) = Fc(:,i,Ny,k,IY)
      end do       ; end do

      do j = 0, Ny ; do i = 0, Nx
         FSz(:,0,i,j) = Fc(:,i,j,0,IZ)
         FSz(:,Nz+1,i,j) = Fc(:,i,j,Nz,IZ)
      end do       ; end do
!
!     Dissipative averaging in complementary points
!     ---------------------------------------------
      do k = 0, Nz ; do j = 0, Ny ; do i = 1, Nx
         call dissipativeFlux(e % storage % Q(:,i-1,j,k), &
                              e % storage % Q(:,i,j,k),   &
                              e % geom % ncXi(:,i,j,k),   &
                              e % geom % t1cXi(:,i,j,k),  &
                              e % geom % t2cXi(:,i,j,k),  &
                              e % geom % JfcXi(i,j,k),    &
                              FVx(:,i,j,k))
      end do                ; end do                ; end do

      do k = 0, Nz ; do i = 0, Nx ; do j = 1, Ny
         call dissipativeFlux(e % storage % Q(:,i,j-1,k), &
                              e % storage % Q(:,i,j,k),   &
                              e % geom % ncEta(:,i,j,k),  &
                              e % geom % t1cEta(:,i,j,k), &
                              e % geom % t2cEta(:,i,j,k), &
                              e % geom % JfcEta(i,j,k),   &
                              FVy(:,j,i,k))
      end do                ; end do                ; end do

      do j = 0, Ny ; do i = 0, Nx ; do k = 1, Nz
         call dissipativeFlux(e % storage % Q(:,i,j,k-1),  &
                              e % storage % Q(:,i,j,k),    &
                              e % geom % ncZeta(:,i,j,k),  &
                              e % geom % t1cZeta(:,i,j,k), &
                              e % geom % t2cZeta(:,i,j,k), &
                              e % geom % JfcZeta(i,j,k),   &
                              FVz(:,k,i,j))
      end do                ; end do                ; end do
!
!     Blending
!     --------
      c2 = (self % c2 * (1.0_RP-switch) + self % c1 * switch)**2
      ! c2 = self % c2 * (1.0_RP-switch) + self % c1 * switch

      do k = 0, Nz ; do j = 0, Ny ; do i = 1, Nx

         p = Pressure(e % storage % Q(:,i-1,j,k))
         invRho = 1.0_RP / e % storage % Q(IRHO,i-1,j,k)
         call getEntropyVariables(e % storage % Q(:,i-1,j,k), p, invRho, w1)

         p = Pressure(e % storage % Q(:,i,j,k))
         invRho = 1.0_RP / e % storage % Q(IRHO,i,j,k)
         call getEntropyVariables(e % storage % Q(:,i,j,k), p, invRho, w2)

         b = dot_product(w2-w1, FSx(:,i,j,k)-FVx(:,i,j,k))
         d = sqrt(b**2 + c2)
         d = (d-b) / d
         ! b = b / c2
         ! if (b <= 0.0_RP) then
         !    d = 1.0_RP
         ! else if (b >= 1.0_RP) then
         !    d = 0.0_RP
         ! else
         !    d = 0.5_RP * (1.0_RP + cos(b*PI))
         ! end if

         FSx(:,i,j,k) = FVx(:,i,j,k) + d*(FSx(:,i,j,k)-FVx(:,i,j,k))

      end do       ; end do       ; end do

      do k = 0, Nz ; do i = 0, Nx ; do j = 1, Ny

         p = Pressure(e % storage % Q(:,i,j-1,k))
         invRho = 1.0_RP / e % storage % Q(IRHO,i,j-1,k)
         call getEntropyVariables(e % storage % Q(:,i,j-1,k), p, invRho, w1)

         p = Pressure(e % storage % Q(:,i,j,k))
         invRho = 1.0_RP / e % storage % Q(IRHO,i,j,k)
         call getEntropyVariables(e % storage % Q(:,i,j,k), p, invRho, w2)

         b = dot_product(w2-w1, FSy(:,j,i,k)-FVy(:,j,i,k))
         d = sqrt(b**2 + c2)
         d = (d-b) / d
         ! b = b / c2
         ! if (b <= 0.0_RP) then
         !    d = 1.0_RP
         ! else if (b >= 1.0_RP) then
         !    d = 0.0_RP
         ! else
         !    d = 0.5_RP * (1.0_RP + cos(b*PI))
         ! end if

         FSy(:,j,i,k) = FVy(:,j,i,k) + d*(FSy(:,j,i,k)-FVy(:,j,i,k))

      end do       ; end do       ; end do

      do j = 0, Ny ; do i = 0, Nx ; do k = 1, Nz

         p = Pressure(e % storage % Q(:,i,j,k-1))
         invRho = 1.0_RP / e % storage % Q(IRHO,i,j,k-1)
         call getEntropyVariables(e % storage % Q(:,i,j,k-1), p, invRho, w1)

         p = Pressure(e % storage % Q(:,i,j,k))
         invRho = 1.0_RP / e % storage % Q(IRHO,i,j,k)
         call getEntropyVariables(e % storage % Q(:,i,j,k), p, invRho, w2)

         b = dot_product(w2-w1, FSz(:,k,i,j)-FVz(:,k,i,j))
         d = sqrt(b**2 + c2)
         d = (d-b) / d
         ! b = b / c2
         ! if (b <= 0.0_RP) then
         !    d = 1.0_RP
         ! else if (b >= 1.0_RP) then
         !    d = 0.0_RP
         ! else
         !    d = 0.5_RP * (1.0_RP + cos(b*PI))
         ! end if

         FSz(:,k,i,j) = FVz(:,k,i,j) + d*(FSz(:,k,i,j)-FVz(:,k,i,j))

      end do       ; end do       ; end do
!
!     Volume integral
!     ---------------
      if (isStrong) then
         div = -ScalarStrongIntegrals % TelescopicVolumeDivergence(e, FSx, FSy, FSz, Fv)
      else
         div = -ScalarWeakIntegrals % TelescopicVolumeDivergence(e, FSx, FSy, FSz, Fv)
      end if

      end associate
      end associate

      computed = .true.

   end function SSFV_hyperbolic
!
!///////////////////////////////////////////////////////////////////////////////
!
   subroutine SSFV_describe(self)
!
!     -------
!     Modules
!     -------
      use MPI_Process_Info, only: MPI_Process
!
!     ---------
!     Interface
!     ---------
      implicit none
      class(SC_SSFV_t), intent(in) :: self


      if (.not. MPI_Process % isRoot) return

      write(STD_OUT,"(30X,A,A30,G0.3,A,G0.3)") "->", "SSFV blending parameter: ", &
                                               self % c1, " - ", self % c2

   end subroutine SSFV_describe
!
!///////////////////////////////////////////////////////////////////////////////
!
!************** Roe-Pike *****************!
!       subroutine dissipativeFlux(QLeft, QRight, nHat, t1, t2, Jf, flux)
!          use FluidData,          only: thermodynamics
!          use VariableConversion, only: getPrimitiveVariables, getRoeVariables
!          use PhysicsStorage
!          implicit none
!          real(kind=RP), intent(in)  :: QLeft(1:NCONS)
!          real(kind=RP), intent(in)  :: QRight(1:NCONS)
!          real(kind=RP), intent(in)  :: nHat(1:NDIM), t1(NDIM), t2(NDIM)
!          real(kind=RP), intent(in)  :: Jf
!          real(kind=RP), intent(out) :: flux(1:NCONS)
! !
! !        ---------------
! !        Local variables
! !        ---------------
! !
!          integer        :: i
!          real(kind=RP)  :: QLRot(5), QRRot(5), VL(NPRIM), VR(NPRIM), aL, aR
!          real(kind=RP)  :: dQ(5), lambda(5), K(5,5), V2abs, alpha(5), dLambda
!          real(kind=RP)  :: rho, u, v, w, V2, H, a
!          real(kind=RP)  :: stab(5)

!          associate(gm1 => thermodynamics % gammaMinus1)
! !
! !        ********************
! !        Perform the rotation
! !        ********************
! !
!          QLRot(1) = QLeft(1)  ; QRRot(1) = QRight(1)

!          QLRot(2) = QLeft (2) * nHat(1) + QLeft (3) * nHat(2) + QLeft (4) * nHat(3)
!          QRRot(2) = QRight(2) * nHat(1) + QRight(3) * nHat(2) + QRight(4) * nHat(3)

!          QLRot(3) = QLeft(2)  * t1(1) + QLeft(3)  * t1(2) + QLeft(4)  * t1(3)
!          QRRot(3) = QRight(2) * t1(1) + QRight(3) * t1(2) + QRight(4) * t1(3)

!          QLRot(4) = QLeft(2)  * t2(1) + QLeft(3)  * t2(2) + QLeft(4)  * t2(3)
!          QRRot(4) = QRight(2) * t2(1) + QRight(3) * t2(2) + QRight(4) * t2(3)

!          QLRot(5) = QLeft(5) ; QRRot(5) = QRight(5)
! !
! !        ***************************
! !        Compute primitive variables
! !        ***************************
! !
!          call getPrimitiveVariables(QLRot, VL)
!          call getPrimitiveVariables(QRRot, VR)

!          aL = sqrt(VL(IPA2))  ; aR = sqrt(VR(IPA2))
! !
! !        *********************
! !        Compute Roe variables: [rho, u, v, w, H, a]
! !        *********************
! !
!          call getRoeVariables(QLRot, QRRot, VL, VR, rho, u, v, w, V2, H, a)
! !
! !        Eigenvalues
! !        -----------
!          lambda(1)   = u-a
!          lambda(2:4) = u
!          lambda(5)   = u+a
! !
! !        Eigenvectors
! !        ------------
!          K(:,1) = (/ 1.0_RP, u-a, v, w, H-u*a /)
!          K(:,2) = (/ 1.0_RP, u, v, w, 0.5_RP*V2 /)
!          K(:,3) = (/ 0.0_RP, 0.0_RP, 1.0_RP, 0.0_RP, v /)
!          K(:,4) = (/ 0.0_RP, 0.0_RP, 0.0_RP, 1.0_RP, w /)
!          K(:,5) = (/ 1.0_RP, u+a, v, w, H+u*a /)
! !
! !        Projections
! !        -----------
!          alpha(1) = ((VR(IPP)-VL(IPP)) - rho * a * (VR(IPU)-VL(IPU)))/(2.0_RP * a * a)
!          alpha(2) = (QRight(IRHO)-QLeft(IRHO)) - (VR(IPP)-VL(IPP))/(a*a)
!          alpha(3) = rho * (VR(IPV)-VL(IPV))
!          alpha(4) = rho * (VR(IPW)-VL(IPW))
!          alpha(5) = ((VR(IPP)-VL(IPP)) + rho * a * (VR(IPU)-VL(IPU)))/(2.0_RP * a * a)
! !
! !        **********************
! !        Perform an entropy fix. Here we use Van Leer's modification of Harten's entropy fix, derived
! !        in: A. Harten, "High resolution schemes for hyperbolic conservation laws". To recover the
! !        Harten entropy fix, set dLambda to 0.5
! !        **********************
! !
! !        Wave #1
! !        -------
!          ! dLambda = max((VR(IPU)-aR) - (VL(IPU)-aL), 0.0_RP)
!          ! if ( abs(lambda(1)) .ge. 2.0_RP * dLambda ) then
!          !    lambda(1) = abs(lambda(1))

!          ! else
!          !    lambda(1) = POW2(lambda(1)) / (4.0_RP * dLambda) + dLambda

!          ! end if
! !
! !        Wave #5
! !        -------
!          ! dLambda = max((VR(IPU)+aR) - (VL(IPU)+aL), 0.0_RP)
!          ! if ( abs(lambda(5)) .ge. 2.0_RP * dLambda ) then
!          !    lambda(5) = abs(lambda(5))

!          ! else
!          !    lambda(5) = POW2(lambda(5)) / (4.0_RP * dLambda) + dLambda

!          ! end if
! !
! !        ****************
! !        Compute the flux
! !        ****************
! !
! !        Perform the average using the averaging function
! !        ------------------------------------------------
!          ! call AveragedStates(QLRot, QRRot, VL(IPP), VR(IPP), VL(IPIRHO), VR(IPIRHO), flux)
!          flux = (EulerFlux1D(QLRot, VL(IPP)) + EulerFlux1D(QRRot, VR(IPP))) * 0.5_RP
! !
! !        Compute the Roe stabilization
! !        -----------------------------
!          ! select case (whichAverage)
!          ! case(PIROZZOLI_SPLIT, KENNEDYGRUBER_SPLIT)
! !
! !           ***************************************************************************
! !           Eigenvalue matrix is corrected for PI and KG variants, see Winters et. al.
! !           "A comparative study on polynomial dealiasing and split form discontinuous
! !           Galerkin schemes for under-resolved turbulence computations"
! !           ***************************************************************************
! !
!             ! lambda(1) = lambda(5)
!          ! end select

!          stab = 0.0_RP
!          do i = 1, 5
!             stab = stab + 0.5_RP * alpha(i) * abs(lambda(i)) * K(:,i)
!          end do
! !
! !        Compute the flux: apply the lambda stabilization here.
! !        ----------------
!          flux = flux - stab
! !
! !        ************************************************
! !        Return momentum equations to the cartesian frame
! !        ************************************************
! !
!          flux(2:4) = nHat*flux(2) + t1*flux(3) + t2*flux(4)
!          flux = flux * Jf

!          end associate
!       end subroutine dissipativeFlux
!************** Low dissipation Roe *****************!
!       subroutine dissipativeFlux(QLeft, QRight, nHat, t1, t2, Jf, flux)
! !
! !        ***********************************************************************
! !           This version, presented by Oßwald et. al. [*], is a modification of
! !           Roe-Pike solver, decreasing the velocity jumps intensity for low
! !           Mach numbers. Normal velocities are scaled such that
! !                         du <- z*du
! !           where z tends to zero as the Mach number tends to zero. Precisely:
! !                         z = min(1, max(ML, MR))
! !
! !           These normal velocities are just scaled to compute Roe dissipation's
! !           intensities (alpha coefficients), not the fluxes (as in the original
! !           Low dissipation method by Thornber et. al.)
! !
! !           * Oßwald et. al. L2Roe: a low dissipation version of Roe’s a
! !              pproximate Riemann solver for low Mach numbers
! !        ***********************************************************************
! !
!          use FluidData, only: thermodynamics

!          implicit none
!          real(kind=RP), intent(in)       :: QLeft(1:NCONS)
!          real(kind=RP), intent(in)       :: QRight(1:NCONS)
!          real(kind=RP), intent(in)       :: nHat(1:NDIM), t1(NDIM), t2(NDIM)
!          real(kind=RP), intent(in)       :: Jf
!          real(kind=RP), intent(out)      :: flux(1:NCONS)
! !
! !        ---------------
! !        Local variables
! !        ---------------
! !
!          integer        :: i
!          real(kind=RP)  :: rhoL, rhouL, rhovL, rhowL, rhoeL, pL, rhoHL, rhoV2L, ML
!          real(kind=RP)  :: rhoR, rhouR, rhovR, rhowR, rhoeR, pR, rhoHR, rhoV2R, MR
!          real(kind=RP)  :: uL, vL, wL, uR, vR, wR, aL, aR, dLambda, z, du, dv, dw
!          real(kind=RP)  :: dp
!          real(kind=RP)  :: QLRot(5), QRRot(5)
!          real(kind=RP)  :: sqrtRhoL, sqrtRhoR, invSumSqrtRhoLR
!          real(kind=RP)  :: invSqrtRhoL, invSqrtRhoR, invRhoL, invRhoR
!          real(kind=RP)  :: rho, u, v, w, H, a, dQ(5), lambda(5), K(5,5), V2abs, alpha(5)
!          real(kind=RP)  :: stab(5)

!          associate(gamma => thermodynamics % gamma, gm1 => thermodynamics % gammaMinus1)
! !
! !        Rotate the variables to the face local frame using normal and tangent vectors
! !        -----------------------------------------------------------------------------
!          rhoL = QLeft(1)                  ; rhoR = QRight(1)
!          invRhoL = 1.0_RP/ rhoL           ; invRhoR = 1.0_RP / rhoR
!          sqrtRhoL = sqrt(rhoL)            ; sqrtRhoR = sqrt(rhoR)
!          invSqrtRhoL = 1.0_RP / sqrtRhoL  ; invSqrtRhoR = 1.0_RP / sqrtRhoR
!          invSumSqrtRhoLR = 1.0_RP / (sqrtRhoL + sqrtRhoR)

!          rhouL = QLeft (2) * nHat(1) + QLeft (3) * nHat(2) + QLeft (4) * nHat(3)
!          rhouR = QRight(2) * nHat(1) + QRight(3) * nHat(2) + QRight(4) * nHat(3)

!          rhovL = QLeft(2)  * t1(1) + QLeft(3)  * t1(2) + QLeft(4)  * t1(3)
!          rhovR = QRight(2) * t1(1) + QRight(3) * t1(2) + QRight(4) * t1(3)

!          rhowL = QLeft(2)  * t2(1) + QLeft(3)  * t2(2) + QLeft(4)  * t2(3)
!          rhowR = QRight(2) * t2(1) + QRight(3) * t2(2) + QRight(4) * t2(3)

!          rhoeL = QLeft(5) ; rhoeR = QRight(5)

!          uL = rhouL * invRhoL    ; uR = rhouR * invRhoR
!          vL = rhovL * invRhoL    ; vR = rhovR * invRhoR
!          wL = rhowL * invRhoL    ; wR = rhowR * invRhoR

!          rhoV2L = (POW2(uL) + POW2(vL) + POW2(wL)) * rhoL
!          rhoV2R = (POW2(uR) + POW2(vR) + POW2(wR)) * rhoR
! !
! !        Compute the enthalpy: here defined as rhoH = gogm1 p + 0.5 rho V^2
! !        --------------------
!          rhoHL = gamma*rhoeL - 0.5_RP*gm1*rhoV2L
!          rhoHR = gamma*rhoeR - 0.5_RP*gm1*rhoV2R

!          pL = gm1 * (rhoeL - 0.5_RP * rhoV2L)
!          pR = gm1 * (rhoeR - 0.5_RP * rhoV2R)

!          aL = sqrt(gamma * pL * invRhoL)
!          aR = sqrt(gamma * pR * invRhoR)
! !
! !        Compute Roe - Pike variables
! !        ----------------------------
!          rho = sqrtRhoL * sqrtRhoR
!          u = (invSqrtRhoL * rhouL + invSqrtRhoR * rhouR) * invSumSqrtRhoLR
!          v = (invSqrtRhoL * rhovL + invSqrtRhoR * rhovR) * invSumSqrtRhoLR
!          w = (invSqrtRhoL * rhowL + invSqrtRhoR * rhowR) * invSumSqrtRhoLR
!          H = (invSqrtRhoL * rhoHL + invSqrtRhoR * rhoHR) * invSumSqrtRhoLR
!          V2abs = POW2(u) + POW2(v) + POW2(w)
!          a = sqrt(gm1*(H - 0.5_RP*V2abs))
! !
! !        Eigenvalues
! !        -----------
!          lambda(1)   = u-a
!          lambda(2:4) = u
!          lambda(5)   = u+a
! !
! !        Eigenvectors
! !        ------------
!          K(:,1) = (/ 1.0_RP, u-a, v, w, H-u*a /)
!          K(:,2) = (/ 1.0_RP, u, v, w, 0.5_RP*V2abs /)
!          K(:,3) = (/ 0.0_RP, 0.0_RP, 1.0_RP, 0.0_RP, v /)
!          K(:,4) = (/ 0.0_RP, 0.0_RP, 0.0_RP, 1.0_RP, w /)
!          K(:,5) = (/ 1.0_RP, u+a, v, w, H+u*a /)
! !
! !        Projections
! !        -----------
! !
! !        ----------------------------------------------------------------------------
! !        Low dissipation Roe-Pike Riemann solver: Reduce the dissipation associated
! !        to the jump in normal velocity. See Obwald et. al. L2Roe: a low dissipation
! !        version of Roe’s approximate Riemann solver for low Mach numbers
! !        ----------------------------------------------------------------------------

!          ML = abs(uL) / aL    ; MR = abs(uR) / aR
!          z  = min(1.0_RP, max(ML,MR))

!          du = z * (uR - uL)
!          dv = z * (vR - vL)
!          dw = z * (wR - wL)
!          dp = pR - pL

!          alpha(1) = (dp - rho * a * du)/(2.0_RP * a * a)
!          alpha(2) = (rhoR-rhoL) - dp/(a*a)
!          alpha(3) = rho * dv
!          alpha(4) = rho * dw
!          alpha(5) = (dp + rho * a * du)/(2.0_RP * a * a)
! !
! !        **********************
! !        Perform an entropy fix. Here we use Van Leer's modification of Harten's entropy fix, derived
! !        in: A. Harten, "High resolution schemes for hyperbolic conservation laws". To recover the
! !        Harten entropy fix, set dLambda to 0.5
! !        **********************
! !
! !        Wave #1
! !        -------
!          ! dLambda = max((uR-aR) - (uL-aL), 0.0_RP)
!          ! if ( abs(lambda(1)) .ge. 2.0_RP * dLambda ) then
!          !    lambda(1) = abs(lambda(1))

!          ! else
!          !    lambda(1) = POW2(lambda(1)) / (4.0_RP * dLambda) + dLambda

!          ! end if
! !
! !        Wave #5
! !        -------
!          ! dLambda = max((uR+aR) - (uL+aL), 0.0_RP)
!          ! if ( abs(lambda(5)) .ge. 2.0_RP * dLambda ) then
!          !    lambda(5) = abs(lambda(5))

!          ! else
!          !    lambda(5) = POW2(lambda(5)) / (4.0_RP * dLambda) + dLambda

!          ! end if
! !
! !        ****************
! !        Compute the flux
! !        ****************
! !
! !        Perform the average using the averaging function
! !        ------------------------------------------------
!          QLRot = (/ rhoL, rhouL, rhovL, rhowL, rhoeL /)
!          QRRot = (/ rhoR, rhouR, rhovR, rhowR, rhoeR /)
!          flux  = (EulerFlux1D(QLRot, pL) + EulerFlux1D(QRRot, pR)) * 0.5_RP
! !
! !        Compute the Roe stabilization
! !        -----------------------------
!          stab = 0.0_RP
!          do i = 1, 5
!             stab = stab + 0.5_RP * alpha(i) * abs(lambda(i)) * K(:,i)
!          end do
! !
! !        Compute the flux: apply the lambda stabilization here.
! !        ----------------
!          flux = flux - stab
! !
! !        ************************************************
! !        Return momentum equations to the cartesian frame
! !        ************************************************
! !
!          flux(2:4) = nHat*flux(2) + t1*flux(3) + t2*flux(4)
!          flux = flux * Jf

!          end associate

!       end subroutine dissipativeFlux
!************** Lax-Friedrichs *****************!
   subroutine dissipativeFlux(Q1, Q2, n, t1, t2, Jf, FV)
!
!     -------
!     Modules
!     -------
      use PhysicsStorage,     only: IRHO, IRHOU, IRHOV, IRHOW, IRHOE
      use VariableConversion, only: Pressure
      use FluidData,          only: thermodynamics
!
!     ---------
!     Interface
!     ---------
      implicit none
      real(RP), intent(in)  :: Q1(NCONS)
      real(RP), intent(in)  :: Q2(NCONS)
      real(RP), intent(in)  :: n(NDIM)
      real(RP), intent(in)  :: t1(NDIM)
      real(RP), intent(in)  :: t2(NDIM)
      real(RP), intent(in)  :: Jf
      real(RP), intent(out) :: FV(NCONS)
!
!     ---------------
!     Local variables
!     ---------------
      real(RP) :: Qn1(NCONS)
      real(RP) :: Qn2(NCONS)
      real(RP) :: F1(NCONS)
      real(RP) :: F2(NCONS)
      real(RP) :: invRho1, invRho2
      real(RP) :: p1, p2
      real(RP) :: u1, u2
      real(RP) :: a1, a2
      real(RP) :: lambda
      real(RP) :: F(NCONS, NDIM)

!
!     Rotate to normal reference frame
!     --------------------------------
      Qn1(IRHO)  = Q1(IRHO)
      Qn1(IRHOU) = Q1(IRHOU)*n(IX)  + Q1(IRHOV)*n(IY)  + Q1(IRHOW)*n(IZ)
      Qn1(IRHOV) = Q1(IRHOU)*t1(IX) + Q1(IRHOV)*t1(IY) + Q1(IRHOW)*t1(IZ)
      Qn1(IRHOW) = Q1(IRHOU)*t2(IX) + Q1(IRHOV)*t2(IY) + Q1(IRHOW)*t2(IZ)
      Qn1(IRHOE) = Q1(IRHOE)

      Qn2(IRHO)  = Q2(IRHO)
      Qn2(IRHOU) = Q2(IRHOU)*n(IX)  + Q2(IRHOV)*n(IY)  + Q2(IRHOW)*n(IZ)
      Qn2(IRHOV) = Q2(IRHOU)*t1(IX) + Q2(IRHOV)*t1(IY) + Q2(IRHOW)*t1(IZ)
      Qn2(IRHOW) = Q2(IRHOU)*t2(IX) + Q2(IRHOV)*t2(IY) + Q2(IRHOW)*t2(IZ)
      Qn2(IRHOE) = Q2(IRHOE)
!
!     Intermediate variables
!     ----------------------
      associate(g => thermodynamics % gamma)
      invRho1 = 1.0_RP / Qn1(IRHO)     ; invRho2 = 1.0_RP / Qn2(IRHO)
      u1      = Qn1(IRHOU) * invRho1   ; u2      = Qn2(IRHOU) * invRho2
      p1      = Pressure(Qn1)          ; p2      = Pressure(Qn2)
      a1      = sqrt(g * p1 * invRho1) ; a2      = sqrt(g * p2 * invRho2)
      end associate
!
!     Average
!     -------
      F1 = EulerFlux1D(Qn1, p1)
      F2 = EulerFlux1D(Qn2, p2)
      FV = AVERAGE(F1, F2)
!
!     Dissipation
!     -----------
      lambda = max(abs(u1)+a1, abs(u2)+a2)
      FV = FV - lambda/2.0_RP * (Qn2-Qn1)

!     Projection to physical reference frame and contravariant scaling
!     ----------------------------------------------------------------
      FV(IRHOU:IRHOW) = FV(IRHOU)*n + FV(IRHOV)*t1 + FV(IRHOW)*t2
      FV = FV * Jf

   end subroutine dissipativeFlux
!
!///////////////////////////////////////////////////////////////////////////////
!
   pure function EulerFlux1D(Q, p) result(F)
!
!     Modules
!     -------
      use PhysicsStorage, only: IRHO, IRHOU, IRHOV, IRHOW, IRHOE
!
!     ---------
!     Interface
!     ---------
      implicit none
      real(RP), intent(in) :: Q(NCONS)
      real(RP), intent(in) :: p
      real(RP)             :: F(NCONS)
!
!     ---------------
!     Local variables
!     ---------------
      real(RP) :: u, v, w

      u = Q(IRHOU) / Q(IRHO)
      v = Q(IRHOV) / Q(IRHO)
      w = Q(IRHOW) / Q(IRHO)

      F(IRHO)  = Q(IRHOU)
      F(IRHOU) = Q(IRHOU) * u + p
      F(IRHOV) = Q(IRHOU) * v
      F(IRHOW) = Q(IRHOU) * w
      F(IRHOE) = ( Q(IRHOE) + p ) * u

   end function EulerFlux1D

end module ShockCapturing
