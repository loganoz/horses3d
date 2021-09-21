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
   use HexMeshClass,               only: HexMesh
   use ElementClass,               only: Element
   use LESModels,                  only: Smagorinsky_t
   use SpectralVanishingViscosity, only: SVV, InitializeSVV
   use DGIntegrals,                only: ScalarWeakIntegrals

   use ShockCapturingKeywords
   use SCsensorClass, only: SCsensor_t, Set_SCsensor, Destruct_SCsensor

   implicit none

   public :: Initialize_ShockCapturing
   public :: ShockCapturingDriver

   type SCdriver_t

         logical :: isActive = .false.  !< On/Off flag
         logical :: hasEllipticTerm     !< .true. if the elliptic term is computed
         logical :: hasHyperbolicTerm   !< .true. if the inviscid term is computed

         real(RP), private :: s1  !< Threshold one (from 0 to 1)
         real(RP), private :: s2  !< Threshold two (from 0 to 1), higher than s1

         type(SCsensor_t),             private              :: sensor     !< Sensor to find discontinuities
         class(ArtificialViscosity_t), private, allocatable :: viscosity  !< Artificial viscosity
         class(ArtificialAdvection_t), private, allocatable :: advection  !< Hyperbolic term

      contains

         procedure :: Detect           => SC_detect
         procedure :: ComputeViscosity => SC_viscosity
         procedure :: ComputeAdvection => SC_advection
         procedure :: Describe         => SC_describe

         final :: SC_destruct

   end type SCdriver_t

   type ArtificialViscosity_t

         integer,  private :: updateMethod     !< Method to compute the viscosity
         integer,  private :: flux_type        !< Art. visc. formulation
         real(RP), private :: mu1              !< First viscosity parameter (low)
         real(RP), private :: alpha1           !< Second viscosity parameter (low)
         real(RP), private :: mu2              !< First viscosity parameter (high)
         real(RP), private :: alpha2           !< Second viscosity parameter (high)
         real(RP), private :: mu2alpha         !< Ratio alpha/mu
         logical,  private :: alphaIsPropToMu  !< .true. if alpha/mu is defined

         type(Smagorinsky_t), private :: Smagorinsky  !< For automatic viscosity

      contains

         procedure :: Initialize => AV_Initialize
         procedure :: Compute    => AV_compute
         procedure :: Describe   => AV_describe

   end type ArtificialViscosity_t

   type ArtificialAdvection_t

      contains

         procedure :: Initialize => AA_Initialize
         procedure :: Compute    => AA_compute
         procedure :: Describe   => AA_describe

   end type ArtificialAdvection_t

   type, extends(ArtificialViscosity_t) :: SC_NoSVV_t

      integer                                          :: fluxType
      procedure(Viscous_Int), nopass, pointer, private :: ViscousFlux => null()

   contains

      procedure :: Initialize => NoSVV_initialize
      procedure :: Compute    => NoSVV_viscosity
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
      procedure :: Compute    => SVV_viscosity
      procedure :: Describe   => SVV_describe

   end type SC_SVV_t

   type, extends(ArtificialAdvection_t) :: SC_SSFV_t

      real(RP) :: c

   contains

      procedure :: Initialize => SSFV_initialize
      procedure :: Compute    => SSFV_hyperbolic
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
   subroutine Initialize_ShockCapturing(self, controlVariables, mesh)
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
      class(HexMesh),                intent(in)    :: mesh
!
!     ---------------
!     Local variables
!     ---------------
      real(RP)                      :: lowThr
      real(RP)                      :: highThr
      real(RP)                      :: thr1
      real(RP)                      :: thr2
      character(len=:), allocatable :: method
      character(len=:), allocatable :: sType
      character(len=:), allocatable :: sVar
      integer                       :: sVarID

!
!     Shock-capturing methods
!     -----------------------
      if (controlVariables % containsKey(SC_VISC_METHOD_KEY)) then
         method = controlVariables % stringValueForKey(SC_VISC_METHOD_KEY, LINE_LENGTH)
      else
         method = SC_NO_VAL
      end if
      call toLower(method)

      allocate(self)

      select case (trim(method))
      case (SC_NOSVV_VAL)
         allocate(SC_NoSVV_t :: self % viscosity)
         self % hasEllipticTerm = .true.

      case (SC_SVV_VAL)
         allocate(SC_SVV_t :: self % viscosity)
         self % hasEllipticTerm = .true.

      case (SC_NO_VAL)
         self % hasEllipticTerm = .false.

      case default
         write(STD_OUT,*) 'ERROR. Unavailable shock-capturing viscous method. Options are:'
         write(STD_OUT,*) '   * ', SC_NO_VAL
         write(STD_OUT,*) '   * ', SC_NOSVV_VAL
         write(STD_OUT,*) '   * ', SC_SVV_VAL
         errorMessage(STD_OUT)
         stop

      end select

      if (controlVariables % containsKey(SC_HYP_METHOD_KEY)) then
         method = controlVariables % stringValueForKey(SC_HYP_METHOD_KEY, LINE_LENGTH)
      else
         method = SC_NO_VAL
      end if
      call toLower(method)

      select case (trim(method))
      case (SC_SSFV_VAL)
         allocate (SC_SSFV_t :: self % advection)
         self % hasHyperbolicTerm = .true.

      case (SC_NO_VAL)
         self % hasHyperbolicTerm = .false.

      case default
         write(STD_OUT,*) 'ERROR. Unavailable shock-capturing hyperbolic method. Options are:'
         write(STD_OUT,*) '   * ', SC_NO_VAL
         write(STD_OUT,*) '   * ', SC_SSFV_VAL
         errorMessage(STD_OUT)
         stop

      end select
!
!     Check if shock-capturing is requested
!     -------------------------------------
      if (controlVariables % containsKey(SC_KEY)) then
         self % isActive = controlVariables%logicalValueForKey(SC_KEY)
      else
         self % isActive = .false.
      end if

      if (.not. self % isActive) then
         deallocate(self)
         if (allocated(self % viscosity)) deallocate(self % viscosity)
         if (allocated(self % advection)) deallocate(self % advection)
         return
      end if
!
!     Initialize viscous and hyperbolic terms
!     ---------------------------------------
      if (allocated(self % viscosity)) call self % viscosity % Initialize(controlVariables, mesh)
      if (allocated(self % advection)) call self % advection % Initialize(controlVariables)
!
!     Sensor thresholds
!     --------------
      if (controlVariables % containsKey(SC_LOW_THRES_KEY)) then
         lowThr = controlVariables % doublePrecisionValueForKey(SC_LOW_THRES_KEY)
      else
         write(STD_OUT,*) 'ERROR. Lower threshold of the sensor must be specified.'
         stop
      end if

      if (controlVariables % containsKey(SC_HIGH_THRES_KEY)) then
         highThr = controlVariables % doublePrecisionValueForKey(SC_HIGH_THRES_KEY)
      else
         write(STD_OUT,*) 'ERROR. Higher threshold of the sensor must be specified.'
         stop
      end if

      if (controlVariables % containsKey(SC_THRES_1_KEY)) then
         self % s1 = controlVariables % doublePrecisionValueForKey(SC_THRES_1_KEY)
      else
         self % s1 = lowThr + (highThr-lowThr) / 3.0_RP
      end if

      if (controlVariables % containsKey(SC_THRES_2_KEY)) then
         self % s2 = controlVariables % doublePrecisionValueForKey(SC_THRES_2_KEY)
      else
         self % s2 = lowThr + (highThr-lowThr) * 2.0_RP / 3.0_RP
      end if
!
!     Sensor type
!     -----------
      if (controlVariables % containsKey(SC_SENSOR_KEY)) then
         sType = controlVariables % stringValueForKey(SC_SENSOR_KEY, LINE_LENGTH)
      else
         sType = SC_GRADRHO_VAL
      end if
      call toLower(sType)
!
!     Sensed variable
!     ---------------
      if (controlVariables % containsKey(SC_VARIABLE_KEY)) then
         sVar = controlVariables % stringValueForKey(SC_VARIABLE_KEY, LINE_LENGTH)
         call toLower(sVar)

         select case (trim(sVar))
         case (SC_RHO_VAL);  sVarID = SC_RHO_ID
         case (SC_RHOU_VAL); sVarID = SC_RHOU_ID
         case (SC_RHOV_VAL); sVarID = SC_RHOV_ID
         case (SC_RHOW_VAL); sVarID = SC_RHOW_ID
         case (SC_RHOE_VAL); sVarID = SC_RHOE_ID
         case (SC_U_VAL);    sVarID = SC_U_ID
         case (SC_V_VAL);    sVarID = SC_V_ID
         case (SC_W_VAL);    sVarID = SC_W_ID
         case (SC_P_VAL);    sVarID = SC_P_ID
         case (SC_RHOP_VAL); sVarID = SC_RHOP_ID
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
         sVarID = SC_RHOP_ID

      end if
!
!     Construct sensor
!     ----------------
      call Set_SCsensor(self % sensor, sType, sVarID, lowThr, highThr)

   end subroutine Initialize_ShockCapturing
!
!///////////////////////////////////////////////////////////////////////////////
!
!     Base class
!
!///////////////////////////////////////////////////////////////////////////////
!
   pure subroutine SC_detect(self, mesh, e)
!
!     ---------
!     Interface
!     ---------
      implicit none
      class(SCdriver_t), intent(in)    :: self
      type(HexMesh),     intent(inout) :: mesh
      type(Element),     intent(inout) :: e


      e % storage % sensor = self % sensor % Compute(mesh, e)

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

!
!     Scale the sensed value to the range [0,1]
!     -----------------------------------------
      switch = self % sensor % Rescale(e % storage % sensor)

      call self % viscosity % Compute(mesh, e, switch, self % hasHyperbolicTerm,  SCflux)

   end subroutine SC_viscosity
!
!///////////////////////////////////////////////////////////////////////////////
!
   function SC_advection(self, e, Fv, Qdot) result(computed)
!
!     ---------
!     Interface
!     ---------
      implicit none
      class(SCdriver_t), intent(in)  :: self
      type(Element),     intent(in)  :: e
      real(RP),          intent(in)  :: Fv(1:NCONS,       &
                                           0:e % Nxyz(1), &
                                           0:e % Nxyz(2), &
                                           0:e % Nxyz(3), &
                                           1:NDIM         )
      real(RP),          intent(out) :: Qdot(1:NCONS,       &
                                             0:e % Nxyz(1), &
                                             0:e % Nxyz(2), &
                                             0:e % Nxyz(3)  )
      logical                        :: computed
!
!     ---------------
!     Local variables
!     ---------------
      real(RP) :: switch

!
!     Scale the sensed value to the range [0,1]
!     -----------------------------------------
      switch = self % sensor % Rescale(e % storage % sensor)

      if (switch >= 1.0_RP) then
         Qdot = self % advection % Compute(e, Fv)
         computed = .true.
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


      if (.not. MPI_Process % isRoot) return

      write(STD_OUT, "(/)")
      call Subsection_Header("Shock-Capturing")

      write(STD_OUT,"(30X,A,A30)", advance="no") "->", "Sensor type: "
      select case (self % sensor % sens_type)
         case (SC_ZERO_ID);    write(STD_OUT,"(A)") SC_ZERO_VAL
         case (SC_ONE_ID);     write(STD_OUT,"(A)") SC_ONE_VAL
         case (SC_MODAL_ID);   write(STD_OUT,"(A)") SC_MODAL_VAL
         case (SC_GRADRHO_ID); write(STD_OUT,"(A)") SC_GRADRHO_VAL
      end select

      write(STD_OUT,"(30X,A,A30)", advance="no") "->", "Sensed variable: "
      select case (self % sensor % sVar)
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

      write(STD_OUT,"(30X,A,A30,F5.1)") "->", "Minimum value: ", self % sensor % low
      write(STD_OUT,"(30X,A,A30,F5.1)") "->", "Threshold s1: ",  self % s1
      write(STD_OUT,"(30X,A,A30,F5.1)") "->", "Threshold s2: ",  self % s2
      write(STD_OUT,"(30X,A,A30,F5.1)") "->", "Maximum value: ", self % sensor % high

      if (allocated(self % viscosity)) call self % viscosity % Describe()
      if (allocated(self % advection)) call self % advection % Describe()

   end subroutine SC_describe
!
!///////////////////////////////////////////////////////////////////////////////
!
   pure subroutine SC_destruct(self)
!
!     ---------
!     Interface
!     ---------
      implicit none
      type(SCdriver_t), intent(inout) :: self


      self%isActive = .false.

      call Destruct_SCsensor(self % sensor)

      if (allocated(self % viscosity)) deallocate(self % viscosity)
      if (allocated(self % advection)) deallocate(self % advection)

   end subroutine SC_destruct
!
!///////////////////////////////////////////////////////////////////////////////
!
!     Artificial viscosity base class
!
!///////////////////////////////////////////////////////////////////////////////
!
   subroutine AV_initialize(self, controlVariables, mesh)
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
      type(HexMesh),                intent(in)    :: mesh
!
!     ---------------
!     Local variables
!     ---------------
      character(len=:), allocatable :: update

!
!     Viscosity values (mu and alpha)
!     -------------------------------
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

      if (controlVariables % containsKey(SC_ALPHA_MU_KEY)) then

         self % alphaIsPropToMu = .true.
         self % mu2alpha        = controlVariables % doublePrecisionValueForKey(SC_ALPHA_MU_KEY)
         self % alpha1          = self % mu2alpha * self % mu1
         self % alpha2          = self % mu2alpha * self % mu2

      else

         self % alphaIsPropToMu = .false.

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

      end if
!
!     Viscosity update method
!     -----------------------
      if (controlVariables % containsKey(SC_VISC_UPDATE_KEY)) then
         update = controlVariables % StringValueForKey(SC_VISC_UPDATE_KEY, LINE_LENGTH)
         call toLower(update)

         select case (trim(update))
         case (SC_CONST_VAL)
            self % updateMethod = SC_CONST_ID

         case (SC_SENSOR_VAL)
            self % updateMethod = SC_SENSOR_ID

         case (SC_SMAG_VAL)

            self % updateMethod = SC_SMAG_ID
            if (.not. self%alphaIsPropToMu) then
               write(STD_OUT,*) 'ERROR. Alpha must be proportional to mu when using SVV-LES.'
               stop
            end if

            ! TODO: Use the default constructor
            self % Smagorinsky % active = .true.
            self % Smagorinsky % requiresWallDistances = .false.
            self % Smagorinsky % WallModel = 0  ! No wall model
            self % Smagorinsky % CS = self % mu1

         case default
            self % updateMethod = SC_CONST_ID

         end select

      else
         self % updateMethod = SC_CONST_ID

      end if

   end subroutine AV_initialize
!
!///////////////////////////////////////////////////////////////////////////////
!
   subroutine AV_compute(self, mesh, e, switch, limit, SCflux)
!
!     ---------
!     Interface
!     ---------
      implicit none
      class(ArtificialViscosity_t), intent(in)    :: self
      type(HexMesh),                intent(inout) :: mesh
      type(Element),                intent(inout) :: e
      real(RP),                     intent(in)    :: switch
      logical,                      intent(in)    :: limit
      real(RP),                     intent(out)   :: SCflux(1:NCONS,       &
                                                            0:e % Nxyz(1), &
                                                            0:e % Nxyz(2), &
                                                            0:e % Nxyz(3), &
                                                            1:NDIM)

      SCflux = 0.0_RP

   end subroutine AV_compute
!
!///////////////////////////////////////////////////////////////////////////////
!
   subroutine AV_describe(self)
!
!     -------
!     Modules
!     -------
      use MPI_Process_Info, only: MPI_Process
      use Headers,          only: Subsection_Header!
!
!     ---------
!     Interface
!     ---------
      implicit none
      class(ArtificialViscosity_t), intent(in) :: self

   end subroutine AV_describe
!
!///////////////////////////////////////////////////////////////////////////////
!
!     Hyperbolic discretization base class
!
!///////////////////////////////////////////////////////////////////////////////
!
   subroutine AA_initialize(self, controlVariables)
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
      class(ArtificialAdvection_t), intent(inout) :: self
      type(FTValueDictionary),      intent(in)    :: controlVariables

   end subroutine AA_initialize
!
!///////////////////////////////////////////////////////////////////////////////
!
   function AA_compute(self, e, Fv) result(div)
!
!     ---------
!     Interface
!     ---------
      implicit none
      class(ArtificialAdvection_t), intent(in) :: self
      type(Element),                intent(in) :: e
      real(RP),                     intent(in) :: Fv(1:NCONS,       &
                                                     0:e % Nxyz(1), &
                                                     0:e % Nxyz(2), &
                                                     0:e % Nxyz(3), &
                                                     1:NDIM         )
      real(RP)                                 :: div(1:NCONS,       &
                                                      0:e % Nxyz(1), &
                                                      0:e % Nxyz(2), &
                                                      0:e % Nxyz(3)  )

      div = 0.0_RP

   end function AA_compute
!
!///////////////////////////////////////////////////////////////////////////////
!
   subroutine AA_describe(self)
!
!     -------
!     Modules
!     -------
      use MPI_Process_Info, only: MPI_Process
      use Headers,          only: Subsection_Header!
!
!     ---------
!     Interface
!     ---------
      implicit none
      class(ArtificialAdvection_t), intent(in) :: self


      if (.not. MPI_Process % isRoot) return

   end subroutine AA_describe
!
!///////////////////////////////////////////////////////////////////////////////
!
!     Simple artificial viscosity
!
!///////////////////////////////////////////////////////////////////////////////
!
   subroutine NoSVV_initialize(self, controlVariables, mesh)
!
!     -------
!     Modules
!     -------
      use FTValueDictionaryClass
      use PhysicsStorage_NS, only: grad_vars, GRADVARS_STATE, &
                                   GRADVARS_ENTROPY, GRADVARS_ENERGY
      use Physics_NS,        only: ViscousFlux_STATE, ViscousFlux_ENTROPY, &
                                   ViscousFlux_ENERGY, GuermondPopovFlux_ENTROPY
!
!     ---------
!     Interface
!     ---------
      implicit none
      class(SC_NoSVV_t),       intent(inout) :: self
      type(FTValueDictionary), intent(in)    :: controlVariables
      type(HexMesh),           intent(in)    :: mesh
!
!     ---------------
!     Local variables
!     ---------------
      character(len=:), allocatable :: flux

!
!     Parent initializer
!     ------------------
      call self % ArtificialViscosity_t % Initialize(controlVariables, mesh)
!
!     Set the flux type
!     -----------------
      if (controlVariables % containsKey(SC_VISC_FLUX_KEY)) then
         flux = controlVariables%stringValueForKey(SC_VISC_FLUX_KEY, LINE_LENGTH)
         call toLower(flux)

         select case (trim(flux))
         case (SC_PHYS_VAL); self % fluxType = SC_PHYS_ID
         case (SC_GP_VAL);   self % fluxType = SC_GP_ID
         case default
            write(STD_OUT,*) 'ERROR. Artificial viscosity type not recognized. Options are:'
            write(STD_OUT,*) '   * ', SC_PHYS_VAL
            write(STD_OUT,*) '   * ', SC_GP_VAL
            errorMessage(STD_OUT)
            stop
         end select

      else
         self % fluxType = SC_PHYS_ID

      end if

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
            errorMessage(STD_OUT)
            stop
         end select

      end select

   end subroutine NoSVV_initialize
!
!///////////////////////////////////////////////////////////////////////////////
!
   subroutine NoSVV_viscosity(self, mesh, e, switch, limit, SCflux)
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
      logical,           intent(in)    :: limit
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


      if ((.not. limit .and. switch > 0.0_RP) .or. &
          (limit .and. switch > 0.0_RP .and. switch < 1.0_RP)) then
!
!        Compute viscosity
!        -----------------
         select case (self % updateMethod)
         case (SC_CONST_ID)
            mu = merge(self % mu2, self % mu1, switch >= 1.0_RP)

         case (SC_SENSOR_ID)
            mu = self % mu1 * (1.0_RP-switch) + self % mu2 * switch

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
         write(STD_OUT,"(30X,A,A30,F4.2)") "->", "LES intensity (CS): ", self % Smagorinsky % CS
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
   subroutine SVV_initialize(self, controlVariables, mesh)
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
      type(HexMesh),           intent(in)    :: mesh

!
!     Parent initializer
!     ------------------
      call self % ArtificialViscosity_t % Initialize(controlVariables, mesh)
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
   subroutine SVV_viscosity(self, mesh, e, switch, limit, SCflux)
!
!     ---------
!     Interface
!     ---------
      implicit none
      class(SC_SVV_t), intent(in)    :: self
      type(HexMesh),   intent(inout) :: mesh
      type(Element),   intent(inout) :: e
      real(RP),        intent(in)    :: switch
      logical,         intent(in)    :: limit
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


      if ((.not. limit .and. switch > 0.0_RP) .or. &
          (limit .and. switch > 0.0_RP .and. switch < 1.0_RP)) then
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

         end select

         if (self % alphaIsPropToMu) then
            sqrt_alpha = self % sqrt_mu2alpha * sqrt_mu
         else
            sqrt_alpha = salpha
         end if
!
!        Compute the viscous flux
!        ------------------------
         call SVV % ComputeInnerFluxes(e, sqrt_mu, sqrt_alpha, switch, SCflux)

      else

         e % storage % SVV_diss = 0.0_RP
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
         write(STD_OUT,"(30X,A,A30,F4.2)") "->", "LES intensity (CS): ", self % Smagorinsky % CS
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
   subroutine SSFV_initialize(self, controlVariables)
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
      class(SC_SSFV_t),        intent(inout) :: self
      type(FTValueDictionary), intent(in)    :: controlVariables


      if (controlVariables % containsKey(SC_BLENDING_KEY)) then
         self % c = controlVariables % doublePrecisionValueForKey(SC_BLENDING_KEY)
      else
         write(STD_OUT,*) "ERROR. A value for the blending parameter must be given."
         stop
      end if

   end subroutine SSFV_initialize
!
!///////////////////////////////////////////////////////////////////////////////
!
   function SSFV_hyperbolic(self, e, Fv) result(div)
!
!     -------
!     Modules
!     -------
      use PhysicsStorage,     only: IRHO
      use NodalStorageClass,  only: NodalStorage
      use Physics,            only: EulerFlux
      use VariableConversion, only: Pressure, getEntropyVariables
!
!     ---------
!     Interface
!     ---------
      implicit none
      class(SC_SSFV_t), intent(in) :: self
      type(Element),    intent(in) :: e
      real(RP),         intent(in) :: Fv(1:NCONS,     &
                                         0:e%Nxyz(1), &
                                         0:e%Nxyz(2), &
                                         0:e%Nxyz(3), &
                                         1:NDIM       )
      real(RP)                     :: div(1:NCONS,     &
                                          0:e%Nxyz(1), &
                                          0:e%Nxyz(2), &
                                          0:e%Nxyz(3)  )
!
!     ---------------
!     Local variables
!     ---------------
      integer  :: i, j, k, r, s
      real(RP) :: FSx(NCONS, 0:e%Nxyz(1)+1, 0:e%Nxyz(2), 0:e%Nxyz(3))
      real(RP) :: FSy(NCONS, 0:e%Nxyz(2)+1, 0:e%Nxyz(1), 0:e%Nxyz(3))
      real(RP) :: FSz(NCONS, 0:e%Nxyz(3)+1, 0:e%Nxyz(1), 0:e%Nxyz(2))
      real(RP) :: FVx(NCONS, 0:e%Nxyz(1)+1, 0:e%Nxyz(2), 0:e%Nxyz(3))
      real(RP) :: FVy(NCONS, 0:e%Nxyz(2)+1, 0:e%Nxyz(1), 0:e%Nxyz(3))
      real(RP) :: FVz(NCONS, 0:e%Nxyz(3)+1, 0:e%Nxyz(1), 0:e%Nxyz(2))
      real(RP) :: F(NCONS, NDIM)
      real(RP) :: w1(NCONS)
      real(RP) :: w2(NCONS)
      real(RP) :: p
      real(RP) :: invRho
      real(RP) :: b
      real(RP) :: d


      associate(Nx => e % Nxyz(1), &
                Ny => e % Nxyz(2), &
                Nz => e % Nxyz(3)  )
      associate(spAxi   => NodalStorage(e % Nxyz(1)), &
                spAeta  => NodalStorage(Ny), &
                spAzeta => NodalStorage(Nz)  )
!
!     Entropy-conservative averaging in complementary points
!     ------------------------------------------------------
      do k = 0, Nz ; do j = 0, Ny ; do i = 1, Nx
         FSx(:,i,j,k) = 0.0_RP
         do s = i, Nx ; do r = 0, i-1
            FSx(:,i,j,k) = FSx(:,i,j,k) + spAxi % Q(r,s) * TwoPointFlux(e % storage % Q(:,r,j,k),    &
                                                                        e % storage % Q(:,s,j,k),    &
                                                                        e % geom % jGradXi(:,r,j,k), &
                                                                        e % geom % jGradXi(:,s,j,k)  )
         end do                  ; end do
         FSx(:,i,j,k) = 2.0_RP * FSx(:,i,j,k)
      end do                ; end do                ; end do

      do k = 0, Nz ; do i = 0, Nx ; do j = 1, Ny
         FSy(:,j,i,k) = 0.0_RP
         do s = j, Ny ; do r = 0, j-1
            FSy(:,j,i,k) = FSy(:,j,i,k) + spAeta % Q(r,s) * TwoPointFlux(e % storage % Q(:,i,r,k),     &
                                                                         e % storage % Q(:,i,s,k),     &
                                                                         e % geom % jGradEta(:,i,r,k), &
                                                                         e % geom % jGradEta(:,i,s,k)  )
         end do                  ; end do
         FSy(:,j,i,k) = 2.0_RP * FSy(:,j,i,k)
      end do                ; end do                ; end do

      do j = 0, Ny ; do i = 1, Nx ; do k = 1, Nz
         FSz(:,k,i,j) = 0.0_RP
         do s = k, Nz ; do r = 0, k-1
            FSz(:,k,i,j) = FSz(:,k,i,j) + spAzeta % Q(r,s) * TwoPointFlux(e % storage % Q(:,i,j,r),      &
                                                                          e % storage % Q(:,i,j,s),      &
                                                                          e % geom % jGradZeta(:,i,j,r), &
                                                                          e % geom % jGradZeta(:,i,j,s)  )
         end do                  ; end do
         FSz(:,k,i,j) = 2.0_RP * FSz(:,k,i,j)
      end do                ; end do                ; end do
!
!     Dissipative averaging in complementary points
!     ---------------------------------------------
      do k = 0, Nz ; do j = 0, Ny ; do i = 1, Nx
         F = dissipativeFlux(e % storage % Q(:,i-1,j,k), e % storage % Q(:,i,j,k))
         FVx(:,i,j,k) = contravariant(F(:,IX), F(:,IY), F(:,IZ), &
                                      e % geom % jGradXi(:,i-1,j,k), e % geom % jGradXi(:,i,j,k))
      end do                ; end do                ; end do

      do k = 0, Nz ; do i = 0, Nx ; do j = 1, Ny
         F = dissipativeFlux(e % storage % Q(:,i,j-1,k), e % storage % Q(:,i,j,k))
         FVy(:,j,i,k) = contravariant(F(:,IX), F(:,IY), F(:,IZ), &
                                      e % geom % jGradEta(:,i,j-1,k), e % geom % jGradEta(:,i,j,k))
      end do                ; end do                ; end do

      do j = 0, Ny ; do i = 1, Nx ; do k = 1, Nz
         F = dissipativeFlux(e % storage % Q(:,i,j,k-1), e % storage % Q(:,i,j,k))
         FVz(:,k,i,j) = contravariant(F(:,IX), F(:,IY), F(:,IZ), &
                                      e % geom % jGradZeta(:,i,j,k-1), e % geom % jGradZeta(:,i,j,k))
      end do                ; end do                ; end do
!
!     Boundaries
!     ----------
      do k = 0, Nz ; do j = 0, Ny

         call EulerFlux(e % storage % Q(:,0,j,k), F)
         FSx(:,0,j,k) = contravariant(F(:,IX), F(:,IY), F(:,IZ), &
                                      e % geom % jGradXi(:,0,j,k), e % geom % jGradXi(:,0,j,k))

         call EulerFlux(e % storage % Q(:,Nx,j,k), F)
         FSx(:,Nx+1,j,k) = contravariant(F(:,IX), F(:,IY), F(:,IZ), &
                                         e % geom % jGradXi(:,Nx,j,k), e % geom % jGradXi(:,Nx,j,k))

      end do       ; end do

      do k = 0, Nz ; do i = 0, Nx

         call EulerFlux(e % storage % Q(:,i,0,k), F)
         FSy(:,0,i,k) = contravariant(F(:,IX), F(:,IY), F(:,IZ), &
                                      e % geom % jGradEta(:,i,0,k), e % geom % jGradEta(:,i,0,k))

         call EulerFlux(e % storage % Q(:,i,Ny,k), F)
         FSy(:,Ny+1,i,k) = contravariant(F(:,IX), F(:,IY), F(:,IZ), &
                                         e % geom % jGradEta(:,i,Ny,k), e % geom % jGradEta(:,i,Ny,k))

      end do       ; end do

      do j = 0, Ny ; do i = 0, Nx

         call EulerFlux(e % storage % Q(:,i,j,0), F)
         FSz(:,0,i,j) = contravariant(F(:,IX), F(:,IY), F(:,IZ), &
                                      e % geom % jGradZeta(:,i,j,0), e % geom % jGradZeta(:,i,j,0))

         call EulerFlux(e % storage % Q(:,i,j,Nz), F)
         FSz(:,Nz+1,i,j) = contravariant(F(:,IX), F(:,IY), F(:,IZ), &
                                         e % geom % jGradZeta(:,i,j,Nz), e % geom % jGradZeta(:,i,j,Nz))

      end do       ; end do
!
!     Blending
!     --------
      do k = 0, Nz ; do j = 0, Ny ; do i = 1, Nx

         p = Pressure(e % storage % Q(:,i-1,j,k))
         invRho = 1.0_RP / e % storage % Q(IRHO,i-1,j,k)
         call getEntropyVariables(e % storage % Q(:,i-1,j,k), p, invRho, w1)

         p = Pressure(e % storage % Q(:,i,j,k))
         invRho = 1.0_RP / e % storage % Q(IRHO,i,j,k)
         call getEntropyVariables(e % storage % Q(:,i,j,k), p, invRho, w2)

         b = dot_product(w2-w1, FSx(:,i,j,k)-FVx(:,i,j,k))
         d = sqrt(b**2 + self % c**2)
         d = (d-b) / d

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
         d = sqrt(b**2 + self % c**2)
         d = (d-b) / d

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
         d = sqrt(b**2 + self % c**2)
         d = (d-b) / d

         FSz(:,k,i,j) = FVz(:,k,i,j) + d*(FSz(:,k,i,j)-FVz(:,k,i,j))

      end do       ; end do       ; end do
!
!     Volume integral
!     ---------------
      div = -ScalarWeakIntegrals % TelescopicVolumeDivergence(e, FSx, FSy, FSz, Fv)

      end associate
      end associate

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

      write(STD_OUT,"(30X,A,A30,G10.3)") "->", "SSFV blending parameter: ", self % c

   end subroutine SSFV_describe
!
!///////////////////////////////////////////////////////////////////////////////
!
   function TwoPointFlux(Qleft, Qright, JaL, JaR) result(F)
!
!     -------
!     Modules
!     -------
      use PhysicsStorage,     only: IRHO, IRHOU, IRHOV, IRHOW
      use VariableConversion, only: Pressure
      use RiemannSolvers_NS,  only: AveragedStates
!
!     ---------
!     Interface
!     ---------
      implicit none
      real(RP), intent(in) :: Qleft(NCONS)
      real(RP), intent(in) :: Qright(NCONS)
      real(RP), intent(in) :: JaL(NDIM)
      real(RP), intent(in) :: JaR(NDIM)
      real(RP)             :: F(NCONS)
!
!     ---------------
!     Local variables
!     ---------------
      real(RP) :: tmp
      real(RP) :: pl, pr
      real(RP) :: iRhoL, iRhoR
      real(RP) :: Ql(NCONS)
      real(RP) :: Qr(NCONS)
      real(RP) :: Fx(NCONS)
      real(RP) :: Fy(NCONS)
      real(RP) :: Fz(NCONS)
      real(RP) :: Ja(NDIM)


      iRhoL = 1.0_RP / Qleft(IRHO) ; iRhoR = 1.0_RP / Qright(IRHO)
      pL    = Pressure(Qleft)      ; pR    = Pressure(Qright)

      Ql = Qleft
      Qr = Qright
      call AveragedStates(Ql, Qr, pL, pR, iRhoL, iRhoR, Fx)

      Ql(IRHOU) = Qleft(IRHOV)  ; Ql(IRHOV) = Qleft(IRHOU)
      Qr(IRHOU) = Qright(IRHOV) ; Qr(IRHOV) = Qright(IRHOU)
      call AveragedStates(Ql, Qr, pL, pR, iRhoL, iRhoR, Fy)
      tmp = Fy(IRHOU) ; Fy(IRHOU) = Fy(IRHOV) ; Fy(IRHOV) = tmp

      Ql = Qleft  ; Ql(IRHOU) = Qleft(IRHOW)  ; Ql(IRHOW) = Qleft(IRHOU)
      Qr = Qright ; Qr(IRHOU) = Qright(IRHOW) ; Qr(IRHOW) = Qright(IRHOU)
      call AveragedStates(Ql, Qr, pL, pR, iRhoL, iRhoR, Fz)
      tmp = Fz(IRHOU) ; Fz(IRHOU) = Fz(IRHOW) ; Fz(IRHOW) = tmp

      F = contravariant(Fx, Fy, Fz, JaL, JaR)

   end function TwoPointFlux
!
!///////////////////////////////////////////////////////////////////////////////
!
   function dissipativeFlux(Q1, Q2) result(FV)
!
!     -------
!     Modules
!     -------
      use Physics,            only: EulerFlux
      use PhysicsStorage,     only: IRHO, IRHOU, IRHOV, IRHOW
      use VariableConversion, only: Pressure
      use FluidData,          only: thermodynamics
!
!     ---------
!     Interface
!     ---------
      implicit none
      real(RP), intent(in) :: Q1(NCONS)
      real(RP), intent(in) :: Q2(NCONS)
      real(RP)             :: FV(NCONS, NDIM)
!
!     ---------------
!     Local variables
!     ---------------
      real(RP) :: F1(NCONS, NDIM)
      real(RP) :: F2(NCONS, NDIM)
      real(RP) :: invRho1, invRho2
      real(RP) :: p1, p2
      real(RP) :: u1, u2
      real(RP) :: v1, v2
      real(RP) :: w1, w2
      real(RP) :: a1, a2
      real(RP) :: lambda
      real(RP) :: jump(NCONS)


      call EulerFlux(Q1, F1)
      call EulerFlux(Q2, F2)
      FV = AVERAGE(F1, F2)

      invRho1 = 1.0_RP / Q1(IRHO)   ; invRho2 = 1.0_RP / Q2(IRHO)
      u1      = Q1(IRHOU) * invRho1 ; u2      = Q2(iRHOU) * invRho2
      v1      = Q1(IRHOV) * invRho1 ; v2      = Q2(iRHOV) * invRho2
      w1      = Q1(IRHOW) * invRho1 ; w2      = Q2(iRHOW) * invRho2
      p1      = Pressure(Q1)        ; p2      = Pressure(Q2)

      associate(g => thermodynamics % gamma)

      jump = Q2 - Q1
      a1 = sqrt(g * p1 * invRho1) ; a2 = sqrt(g * p2 * invRho2)

      lambda = max(abs(u1)+a1, abs(u2)+a2)
      FV(:,IX) = FV(:,IX) - lambda/2.0_RP * jump

      lambda = max(abs(v1)+a1, abs(v2)+a2)
      FV(:,IY) = FV(:,IY) - lambda/2.0_RP * jump

      lambda = max(abs(w1)+a1, abs(w2)+a2)
      FV(:,IZ) = FV(:,IZ) - lambda/2.0_RP * jump

      end associate

   end function dissipativeFlux
!
!///////////////////////////////////////////////////////////////////////////////
!
   pure function contravariant(Fx, Fy, Fz, JaL, JaR) result(Fcont)
!
!     ---------
!     Interface
!     ---------
      implicit none
      real(RP), intent(in) :: Fx(NCONS)
      real(RP), intent(in) :: Fy(NCONS)
      real(RP), intent(in) :: Fz(NCONS)
      real(RP), intent(in) :: JaL(NDIM)
      real(RP), intent(in) :: JaR(NDIM)
      real(RP)             :: Fcont(NCONS)
!
!     ---------------
!     Local variables
!     ---------------
      real(RP) :: Ja(NDIM)


      Ja = AVERAGE(JaL, JaR)
      Fcont = Fx*Ja(IX) + Fy*Ja(IY) + Fz*Ja(IZ)

   end function contravariant

end module ShockCapturing
