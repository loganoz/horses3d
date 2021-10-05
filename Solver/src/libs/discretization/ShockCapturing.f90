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
   use LESModels,                  only: Smagorinsky_t
   use SpectralVanishingViscosity, only: SVV, InitializeSVV
   use DGIntegrals,                only: ScalarWeakIntegrals, ScalarStrongIntegrals

   use ShockCapturingKeywords
   use SCsensorClass, only: SCsensor_t, Set_SCsensor, Destruct_SCsensor

   implicit none

   public :: Initialize_ShockCapturing
   public :: ShockCapturingDriver

   type SCdriver_t

         logical :: isActive = .false.  !< On/Off flag
         logical :: hasEllipticTerm     !< .true. if the elliptic term is computed
         logical :: hasHyperbolicTerm   !< .true. if the inviscid term is computed

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

      real(RP) :: c1  ! Lower limit of the blending parameter ("More FV")
      real(RP) :: c2  ! Higher limit of the blending parameter ("More FS")

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
         updateCompGeometry = .true.

      case (SC_NO_VAL)
         self % hasHyperbolicTerm = .false.
         updateCompGeometry = .false.

      case default
         write(STD_OUT,*) 'ERROR. Unavailable shock-capturing hyperbolic method. Options are:'
         write(STD_OUT,*) '   * ', SC_NO_VAL
         write(STD_OUT,*) '   * ', SC_SSFV_VAL
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
         return
      end if
!
!     Initialize viscous and hyperbolic terms
!     ---------------------------------------
      if (allocated(self % viscosity)) call self % viscosity % Initialize(controlVariables, sem % mesh)
      if (allocated(self % advection)) call self % advection % Initialize(controlVariables, sem % mesh)
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

!
!     Scale the sensed value to the range [0,1]
!     -----------------------------------------
      switch = self % sensor % Rescale(e % storage % sensor)

      call self % viscosity % Compute(mesh, e, switch, SCflux)

   end subroutine SC_viscosity
!
!///////////////////////////////////////////////////////////////////////////////
!
   subroutine SC_advection(self, e, Fv, Qdot, isStrong)
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
!
!     ---------------
!     Local variables
!     ---------------
      real(RP) :: switch
      logical  :: strongIntegral

!
!     Scale the sensed value to the range [0,1]
!     -----------------------------------------
      switch = self % sensor % Rescale(e % storage % sensor)

      if (present(isStrong)) then
         strongIntegral = isStrong
      else
         strongIntegral = .false.
      end if

      Qdot = self % advection % Compute(e, switch, Fv, strongIntegral)

   end subroutine SC_advection
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

      call self % sensor % Describe()
      if (allocated(self % viscosity)) call self % viscosity % Describe()
      if (allocated(self % advection)) call self % advection % Describe()

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
   subroutine AV_compute(self, mesh, e, switch, SCflux)
!
!     ---------
!     Interface
!     ---------
      implicit none
      class(ArtificialViscosity_t), intent(in)    :: self
      type(HexMesh),                intent(inout) :: mesh
      type(Element),                intent(inout) :: e
      real(RP),                     intent(in)    :: switch
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
   subroutine AA_initialize(self, controlVariables, mesh)
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
      type(HexMesh),                intent(inout) :: mesh

   end subroutine AA_initialize
!
!///////////////////////////////////////////////////////////////////////////////
!
   function AA_compute(self, e, switch, Fv, isStrong) result(div)
!
!     ---------
!     Interface
!     ---------
      implicit none
      class(ArtificialAdvection_t), intent(in) :: self
      type(Element),                intent(in) :: e
      real(RP),                     intent(in) :: switch
      logical,                      intent(in) :: isStrong
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
   subroutine SSFV_initialize(self, controlVariables, mesh)
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
      class(SC_SSFV_t),        intent(inout) :: self
      type(FTValueDictionary), intent(in)    :: controlVariables
      type(HexMesh),           intent(inout) :: mesh
!
!     ---------------
!     Local variables
!     ---------------
!
      integer                          :: eID
      type(TransfiniteHexMap), pointer :: hexMap    => null()
      type(TransfiniteHexMap), pointer :: hex8Map   => null()
      type(TransfiniteHexMap), pointer :: genHexMap => null()

!
!     Blending parameters
!     -------------------
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
   function SSFV_hyperbolic(self, e, switch, Fv, isStrong) result(div)
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
      real(RP),         intent(in) :: switch
      logical,          intent(in) :: isStrong
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
      real(RP) :: c
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
!     Boundaries
!     ----------
      do k = 0, Nz ; do j = 0, Ny

         call EulerFlux(e % storage % Q(:,0,j,k), F)
         FSx(:,0,j,k) = contravariant(F(:,IX), F(:,IY), F(:,IZ),   &
                                      e % geom % jGradXi(:,0,j,k), &
                                      e % geom % jGradXi(:,0,j,k))

         call EulerFlux(e % storage % Q(:,Nx,j,k), F)
         FSx(:,Nx+1,j,k) = contravariant(F(:,IX), F(:,IY), F(:,IZ),    &
                                         e % geom % jGradXi(:,Nx,j,k), &
                                         e % geom % jGradXi(:,Nx,j,k))

      end do       ; end do

      do k = 0, Nz ; do i = 0, Nx

         call EulerFlux(e % storage % Q(:,i,0,k), F)
         FSy(:,0,i,k) = contravariant(F(:,IX), F(:,IY), F(:,IZ),    &
                                      e % geom % jGradEta(:,i,0,k), &
                                      e % geom % jGradEta(:,i,0,k))

         call EulerFlux(e % storage % Q(:,i,Ny,k), F)
         FSy(:,Ny+1,i,k) = contravariant(F(:,IX), F(:,IY), F(:,IZ),     &
                                         e % geom % jGradEta(:,i,Ny,k), &
                                         e % geom % jGradEta(:,i,Ny,k))

      end do       ; end do

      do j = 0, Ny ; do i = 0, Nx

         call EulerFlux(e % storage % Q(:,i,j,0), F)
         FSz(:,0,i,j) = contravariant(F(:,IX), F(:,IY), F(:,IZ),     &
                                      e % geom % jGradZeta(:,i,j,0), &
                                      e % geom % jGradZeta(:,i,j,0))

         call EulerFlux(e % storage % Q(:,i,j,Nz), F)
         FSz(:,Nz+1,i,j) = contravariant(F(:,IX), F(:,IY), F(:,IZ),      &
                                         e % geom % jGradZeta(:,i,j,Nz), &
                                         e % geom % jGradZeta(:,i,j,Nz))

      end do       ; end do
!
!     Add dissipation if required only
!     --------------------------------
      if (switch > 0.0_RP) then
!
!        Dissipative averaging in complementary points
!        ---------------------------------------------
         do k = 0, Nz ; do j = 0, Ny ; do i = 1, Nx
            FVx(:,i,j,k) = dissipativeFlux(e % storage % Q(:,i-1,j,k), &
                                          e % storage % Q(:,i,j,k),   &
                                          e % geom % ncXi(:,i,j,k),   &
                                          e % geom % t1cXi(:,i,j,k),  &
                                          e % geom % t2cXi(:,i,j,k),  &
                                          e % geom % JfcXi(i,j,k))
         end do                ; end do                ; end do

         do k = 0, Nz ; do i = 0, Nx ; do j = 1, Ny
            FVy(:,j,i,k) = dissipativeFlux(e % storage % Q(:,i,j-1,k), &
                                          e % storage % Q(:,i,j,k),   &
                                          e % geom % ncEta(:,i,j,k),  &
                                          e % geom % t1cEta(:,i,j,k), &
                                          e % geom % t2cEta(:,i,j,k), &
                                          e % geom % JfcEta(i,j,k))
         end do                ; end do                ; end do

         do j = 0, Ny ; do i = 1, Nx ; do k = 1, Nz
            FVz(:,k,i,j) = dissipativeFlux(e % storage % Q(:,i,j,k-1),  &
                                          e % storage % Q(:,i,j,k),    &
                                          e % geom % ncZeta(:,i,j,k),  &
                                          e % geom % t1cZeta(:,i,j,k), &
                                          e % geom % t2cZeta(:,i,j,k), &
                                          e % geom % JfcZeta(i,j,k))
         end do                ; end do                ; end do
!
!        Blending
!        --------
         c = self % c2 * (1.0_RP-switch) + self % c1 * switch

         do k = 0, Nz ; do j = 0, Ny ; do i = 1, Nx

            p = Pressure(e % storage % Q(:,i-1,j,k))
            invRho = 1.0_RP / e % storage % Q(IRHO,i-1,j,k)
            call getEntropyVariables(e % storage % Q(:,i-1,j,k), p, invRho, w1)

            p = Pressure(e % storage % Q(:,i,j,k))
            invRho = 1.0_RP / e % storage % Q(IRHO,i,j,k)
            call getEntropyVariables(e % storage % Q(:,i,j,k), p, invRho, w2)

            b = dot_product(w2-w1, FSx(:,i,j,k)-FVx(:,i,j,k))
            d = sqrt(b**2 + c**2)
            d = (d-b) / d
            d = min(d, 1.0_RP)

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
            d = sqrt(b**2 + c**2)
            d = (d-b) / d
            d = min(d, 1.0_RP)

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
            d = sqrt(b**2 + c**2)
            d = (d-b) / d
            d = min(d, 1.0_RP)

            FSz(:,k,i,j) = FVz(:,k,i,j) + d*(FSz(:,k,i,j)-FVz(:,k,i,j))

         end do       ; end do       ; end do

      end if
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

      write(STD_OUT,"(30X,A,A30,G10.3,A,G10.3)") "->", "SSFV blending parameter: ", &
                                                 self % c1, " - ", self % c2

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
   function dissipativeFlux(Q1, Q2, n, t1, t2, Jf) result(FV)
!
!     -------
!     Modules
!     -------
      use Physics,            only: EulerFlux
      use PhysicsStorage,     only: IRHO, IRHOU, IRHOV, IRHOW, IRHOE
      use VariableConversion, only: Pressure
      use FluidData,          only: thermodynamics
!
!     ---------
!     Interface
!     ---------
      implicit none
      real(RP), intent(in) :: Q1(NCONS)
      real(RP), intent(in) :: Q2(NCONS)
      real(RP), intent(in) :: n(NDIM)
      real(RP), intent(in) :: t1(NDIM)
      real(RP), intent(in) :: t2(NDIM)
      real(RP), intent(in) :: Jf
      real(RP)             :: FV(NCONS)
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
      u1      = Qn1(IRHOU) * invRho1   ; u2      = Qn2(iRHOU) * invRho2
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
!
!     Projection to physical reference frame and contravariant scaling
!     ----------------------------------------------------------------
      FV(IRHOU:IRHOW) = FV(IRHOU)*n + FV(IRHOV)*t1 + FV(IRHOW)*t2
      FV = FV * Jf

   end function dissipativeFlux
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
!
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
