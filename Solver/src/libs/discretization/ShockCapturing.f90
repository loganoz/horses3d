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
   use PhysicsStorage_NS,          only: NCONS, NGRAD
   use FluidData_NS,               only: dimensionless
   use Utilities,                  only: toLower
   use HexMeshClass,               only: HexMesh
   use ElementClass,               only: Element
   use LESModels,                  only: Smagorinsky_t
   use SpectralVanishingViscosity, only: SVV, InitializeSVV

   use ShockCapturingKeywords
   use SCsensorClass, only: SCsensor_t, Set_SCsensor, Destruct_SCsensor

   implicit none

   public :: Initialize_ShockCapturing
   public :: ShockCapturingDriver

   type SCdriver_t

         logical :: isActive = .false.  !< On/Off flag

         logical          :: hasHyperbolicTerm   !< .true. if the inviscid term is computed
         integer, private :: updateMethod        !< Method to compute the viscosity
         integer, private :: flux_type           !< Art. visc. formulation

         real(RP), private :: s1     !< Threshold one (from 0 to 1)
         real(RP), private :: s2     !< Threshold two (from 0 to 1), higher that s1

         real(RP), private :: mu1              !< First viscosity parameter (low)
         real(RP), private :: alpha1           !< Second viscosity parameter (low)
         real(RP), private :: mu2              !< First viscosity parameter (high)
         real(RP), private :: alpha2           !< Second viscosity parameter (high)
         real(RP), private :: mu2alpha         !< Ratio alpha/mu
         logical,  private :: alphaIsPropToMu  !< .true. if alpha/mu is defined

         type(SCsensor_t), private :: sensor  !< Sensor to find discontinuities

      contains

         procedure :: Initialize => SC_initialize
         procedure :: Detect     => SC_detect
         procedure :: Viscosity  => SC_viscosity
         !procedure :: Hyperbolic => SC_hyperbolic
         procedure :: Describe   => SC_describe

         final :: SC_destruct

   end type SCdriver_t

   type, extends(SCdriver_t) :: ArtViscDriver_t
      integer                                          :: fluxType
      procedure(Viscous_Int), nopass, pointer, private :: ViscousFlux => null()
   contains
      procedure :: Initialize => ArtVisc_initialize
      procedure :: Viscosity  => ArtVisc_viscosity
      procedure :: Describe   => ArtVisc_describe
      final :: ArtVisc_destruct
   end type ArtViscDriver_t

   type, extends(SCdriver_t) :: SVVdriver_t
      real(RP), private :: sqrt_mu1
      real(RP), private :: sqrt_alpha1
      real(RP), private :: sqrt_mu2
      real(RP), private :: sqrt_alpha2
      real(RP), private :: sqrt_mu2alpha
   contains
      procedure :: Initialize => SVV_Initialize
      procedure :: Viscosity  => SVV_viscosity
      procedure :: Describe   => SVV_describe
   end type SVVdriver_t

   type, extends(SVVdriver_t) :: SSFVdriver_t
      real(RP) :: c
   contains
      !procedure :: Initialize => SSFV_initialize
      !procedure :: Viscosity  => SSFV_viscosity
      !procedure :: Hyperbolic => SSFV_hyperbolic
      !procedure :: Describe   => SSFV_describe
   end type SSFVdriver_t
!
!  Interfaces
!  ----------
   abstract interface
      pure subroutine Viscous_Int(nEqn, nGradEqn, Q, Q_x, Q_y, Q_z, mu, beta, kappa, F)
         import RP, NDIM
         integer,       intent(in)  :: nEqn
         integer,       intent(in)  :: nGradEqn
         real(kind=RP), intent(in)  :: Q   (1:nEqn     )
         real(kind=RP), intent(in)  :: Q_x (1:nGradEqn)
         real(kind=RP), intent(in)  :: Q_y (1:nGradEqn)
         real(kind=RP), intent(in)  :: Q_z (1:nGradEqn)
         real(kind=RP), intent(in)  :: mu
         real(kind=RP), intent(in)  :: beta
         real(kind=RP), intent(in)  :: kappa
         real(kind=RP), intent(out) :: F(1:nEqn, 1:NDIM)
      end subroutine Viscous_Int
   end interface

   class(SCdriver_t), allocatable :: ShockCapturingDriver
   type(Smagorinsky_t)            :: Smagorinsky
!
!  ========
   contains
!  ========
!
!
!///////////////////////////////////////////////////////////////////////////////
!
!     Initializers
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
      class(SCdriver_t), allocatable, intent(inout) :: self
      class(FTValueDictionary),       intent(in)    :: controlVariables
      class(HexMesh),                 intent(in)    :: mesh
!
!     ---------------
!     Local variables
!     ---------------
      real(RP)                      :: lowThr
      real(RP)                      :: highThr
      real(RP)                      :: thr1
      real(RP)                      :: thr2
      character(len=:), allocatable :: method
      character(len=:), allocatable :: update
      character(len=:), allocatable :: sType
      character(len=:), allocatable :: sVar
      integer                       :: sVarID

!
!     Shock-capturing method
!     ----------------------
      if (controlVariables % containsKey(SC_METHOD_KEY)) then
         method = controlVariables % stringValueForKey(SC_METHOD_KEY, LINE_LENGTH)
      else
         method = SC_SSFV_VAL
      end if
      call toLower(method)

      select case (trim(method))
      case (SC_NOSVV_VAL)
         allocate(ArtViscDriver_t :: self)
         self % hasHyperbolicTerm = .false.

      case (SC_SVV_VAL)
         allocate(SVVdriver_t :: self)
         self % hasHyperbolicTerm = .false.

      case (SC_SSFV_VAL)
         allocate(SSFVdriver_t :: self)
         self % hasHyperbolicTerm = .true.

      case default
         write(STD_OUT,*) 'ERROR. Unavailable shock-capturing method. Options are:'
         write(STD_OUT,*) '   * ', SC_NOSVV_VAL
         write(STD_OUT,*) '   * ', SC_SVV_VAL
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

      if (.not. self % isActive) return
!
!     Viscosity values (mu and alpha)
!     -------------------------------
      if (controlVariables % containsKey(SC_MU1_KEY)) then
         self % mu1 = controlVariables % doublePrecisionValueForKey(SC_MU1_KEY)
      else
         write(STD_OUT,*) "ERROR. A value for the artificial 'mu 1' must be given."
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
      if (controlVariables % containsKey(SC_UPDATE_KEY)) then
         update = controlVariables % StringValueForKey(SC_UPDATE_KEY, LINE_LENGTH)
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
            Smagorinsky % active = .true.
            Smagorinsky % requiresWallDistances = .false.
            Smagorinsky % WallModel = 0  ! No wall model
            Smagorinsky % CS = self % mu1

         case default
            self % updateMethod = SC_CONST_ID

         end select

      else
         self % updateMethod = SC_CONST_ID

      end if
!
!     Implementation-specific intialization
!     -------------------------------------
      call self % Initialize(controlVariables, mesh)
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
   subroutine SC_initialize(self, controlVariables, mesh)
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
      class(SCdriver_t),       intent(inout) :: self
      type(FTValueDictionary), intent(in)    :: controlVariables
      type(HexMesh),           intent(in)    :: mesh


   end subroutine SC_initialize
!
!///////////////////////////////////////////////////////////////////////////////
!
   subroutine ArtVisc_initialize(self, controlVariables, mesh)
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
      class(ArtViscDriver_t),  intent(inout) :: self
      type(FTValueDictionary), intent(in)    :: controlVariables
      type(HexMesh),           intent(in)    :: mesh
!
!     ---------------
!     Local variables
!     ---------------
      character(len=:), allocatable :: flux


      ! Set the flux type
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

      ! Now point to the correct function
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

   end subroutine ArtVisc_initialize
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
      class(SVVdriver_t),      intent(inout) :: self
      type(FTValueDictionary), intent(in)    :: controlVariables
      type(HexMesh),           intent(in)    :: mesh


      ! Set the square root of the viscosities
      self % sqrt_mu1 = sqrt(self % mu1)
      self % sqrt_mu2 = sqrt(self % mu2)

      if (self % alphaIsPropToMu) then
         self % sqrt_mu2alpha = sqrt(self % mu2alpha)
      else
         self % sqrt_alpha1 = sqrt(self % alpha1)
         self % sqrt_alpha2 = sqrt(self % alpha2)
      end if

      ! Start the SVV module
      call InitializeSVV(SVV, controlVariables, mesh)

   end subroutine SVV_initialize
!
!///////////////////////////////////////////////////////////////////////////////
!
!     Destructors
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

   end subroutine SC_destruct
!
!///////////////////////////////////////////////////////////////////////////////
!
   pure subroutine ArtVisc_destruct(self)
!
!     ---------
!     Interface
!     ---------
      implicit none
      type(ArtViscDriver_t), intent(inout) :: self


      if (associated(self % ViscousFlux)) nullify(self % ViscousFlux)

   end subroutine ArtVisc_destruct
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
      real(RP),          intent(out)   :: &
         SCflux(1:NCONS, 0:e % Nxyz(1), 0:e % Nxyz(2), 0:e % Nxyz(3), 1:NDIM)


      SCflux = 0.0_RP

   end subroutine SC_viscosity
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

      write(STD_OUT,"(30X,A,A30)", advance="no") "->", "Viscosity update method: "
      select case (self % updateMethod)
         case (SC_CONST_ID);  write(STD_OUT,"(A)") SC_CONST_VAL
         case (SC_SENSOR_ID); write(STD_OUT,"(A)") SC_SENSOR_VAL
         case (SC_SMAG_ID);   write(STD_OUT,"(A)") SC_SMAG_VAL
      end select

      if (self % updateMethod == SC_SMAG_ID) then
         write(STD_OUT,"(30X,A,A30,F4.2)") "->", "LES intensity (CS): ", Smagorinsky % CS
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

      ! The rest goes in the derived clases

   end subroutine SC_describe
!
!///////////////////////////////////////////////////////////////////////////////
!
!     Simple artificial viscosity
!
!///////////////////////////////////////////////////////////////////////////////
!
   subroutine ArtVisc_viscosity(self, mesh, e, SCflux)
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
      class(ArtViscDriver_t), intent(in)    :: self
      type(HexMesh),          intent(inout) :: mesh
      type(Element),          intent(inout) :: e
      real(RP),               intent(out)   :: SCflux(1:NCONS,       &
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
      real(RP) :: switch
      real(RP) :: delta
      real(RP) :: kappa
      real(RP) :: mu(0:e % Nxyz(1), 0:e % Nxyz(2), 0:e % Nxyz(3))
      real(RP) :: covariantFlux(1:NCONS, 1:NDIM)

!
!     Scale the sensed value to the range [0,1]
!     -----------------------------------------
      switch = self % sensor % Rescale(e % storage % sensor)

      if (switch > 0.0_RP) then
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
               call Smagorinsky % ComputeViscosity(delta, e % geom % dWall(i,j,k), &
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

   end subroutine ArtVisc_viscosity
!
!///////////////////////////////////////////////////////////////////////////////
!
   subroutine ArtVisc_describe(self)
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
      class(ArtViscDriver_t), intent(in) :: self


      if (.not. MPI_Process % isRoot) return

      call self % SCdriver_t % Describe()

      write(STD_OUT,"(30X,A,A30)", advance="no") "->", "Dissipation type: "
      select case (self % fluxType)
         case (SC_PHYS_ID); write(STD_OUT,"(A)") SC_PHYS_VAL
         case (SC_GP_ID);   write(STD_OUT,"(A)") SC_GP_VAL
      end select

   end subroutine ArtVisc_describe
!
!///////////////////////////////////////////////////////////////////////////////
!
!     SVV filtered viscosity
!
!///////////////////////////////////////////////////////////////////////////////
!
   subroutine SVV_viscosity(self, mesh, e, SCflux)
!
!     ---------
!     Interface
!     ---------
      implicit none
      class(SVVdriver_t), intent(in)    :: self
      type(HexMesh),      intent(inout) :: mesh
      type(Element),      intent(inout) :: e
      real(RP),           intent(out)   :: SCflux(1:NCONS,       &
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
      real(RP) :: switch
      real(RP) :: delta
      real(RP) :: salpha
      real(RP) :: sqrt_mu(0:e % Nxyz(1), 0:e % Nxyz(2), 0:e % Nxyz(3))
      real(RP) :: sqrt_alpha(0:e % Nxyz(1), 0:e % Nxyz(2), 0:e % Nxyz(3))

!
!     Scale the sensed value to the range [0,1]
!     -----------------------------------------
      switch = self % sensor % Rescale(e % storage % sensor)

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

         case (SC_SMAG_ID)

            delta = (e % geom % Volume / product(e % Nxyz + 1)) ** (1.0_RP / 3.0_RP)
            do k = 0, e % Nxyz(3) ; do j = 0, e % Nxyz(2) ; do i = 0, e % Nxyz(1)
               call Smagorinsky % ComputeViscosity(delta, e % geom % dWall(i,j,k), &
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
         call SVV % ComputeInnerFluxes(mesh, e, sqrt_mu, sqrt_alpha, switch, SCflux)

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
      class(SVVdriver_t), intent(in) :: self


      if (.not. MPI_Process % isRoot) return

      call self % SCdriver_t % Describe()
      call SVV % Describe()

   end subroutine SVV_describe

end module ShockCapturing
