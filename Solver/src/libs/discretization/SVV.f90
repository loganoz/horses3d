!
!//////////////////////////////////////////////////////
!
!   @File:    SVV.f90
!   @Author:  Juan Manzanero (juan.manzanero@upm.es)
!   @Created: Sat Jan  6 11:47:48 2018
!   @Last revision date: Tue Aug 31 18:52:22 2021
!   @Last revision author: Wojciech Laskowski (wj.laskowski@upm.es)
!   @Last revision commit: f6acb799fc59cc2a7f6930d9a86a5f7527bc5ced
!
!//////////////////////////////////////////////////////
!
#include "Includes.h"
#if defined(NAVIERSTOKES)
module SpectralVanishingViscosity
   use SMConstants
   use MeshTypes
   use Physics
   use PhysicsStorage
   use MPI_Face_Class
   use EllipticDiscretizationClass
   use HexMeshClass
   use NodalStorageClass
   use GaussQuadrature
   use FluidData
   use MPI_Process_Info             , only: MPI_Process
   use Headers                      , only: Subsection_Header
   use LESModels                    , only: Smagorinsky_t
   use Utilities                    , only: toLower
   implicit none

   private
   public   SVV, InitializeSVV

   integer,          parameter  :: Nmax = 20
!
!  Keywords
!  --------
   character(len=*), parameter  :: SVV_KEY               =  "enable svv"
   character(len=*), parameter  :: SVV_MU_KEY            =  "svv viscosity"
   character(len=*), parameter  :: SVV_ALPHA_KEY         =  "svv alpha viscosity"
   character(len=*), parameter  :: SVV_ALPHA_MU_KEY      =  "svv alpha/mu ratio"
   character(len=*), parameter  :: SVV_CUTOFF_KEY        =  "svv filter cutoff"
   character(len=*), parameter  :: SVV_MU_SMAG           =  "svv les intensity"
   character(len=*), parameter  :: FILTER_SHAPE_KEY      =  "svv filter shape"
   character(len=*), parameter  :: FILTER_TYPE_KEY       =  "svv filter type"
   character(len=*), parameter  :: DISSIPATION_TYPE_KEY  =  "svv dissipation type"
!
!  Filter types
!  ------------
   enum, bind(C)
      enumerator :: HPASS_FILTER, LPASS_FILTER
   end enum
!
!  Filter shapes
!  -------------
   enum, bind(C)
      enumerator :: POW_FILTER, SHARP_FILTER, EXP_FILTER
   end enum
!
!  Dissipation types
!  -----------------
   enum, bind(C)
      enumerator :: PHYSICAL_DISS, GUERMOND_DISS
   end enum

   type FilterMatrices_t
      logical                    :: constructed = .false.
      integer                    :: N
      real(kind=RP), allocatable :: Q(:,:)
      real(kind=RP), allocatable :: F(:,:)
      real(kind=RP), allocatable :: B(:,:)
      contains
         procedure :: Recompute => FilterMatrices_Recompute
   end type FilterMatrices_t

   type  SVV_t
      logical                                     :: automatic_Psvv
      logical                                     :: enabled
      logical                                     :: muIsSmagorinsky = .FALSE.
      logical                                     :: alphaIsPropToMu = .FALSE.
      integer                                     :: filterType
      integer                                     :: filterShape
      integer                                     :: diss_type
      integer, allocatable                        :: entropy_indexes(:)
      real(kind=RP)                               :: muSVV,    sqrt_muSVV
      real(kind=RP)                               :: alphaSVV, sqrt_alphaSVV
      real(kind=RP)                               :: sqrt_mu2alpha
      real(kind=RP)                               :: Psvv
      type(FilterMatrices_t)                      :: filters(0:Nmax)
      procedure(Compute_Hflux_f), nopass, pointer :: Compute_Hflux
      procedure(Compute_SVV_f),   nopass, pointer :: Compute_SVV
      contains
         procedure      :: ConstructFilter    => SVV_ConstructFilter
         procedure      :: ComputeInnerFluxes => SVV_ComputeInnerFluxes
         procedure      :: Describe           => SVV_Describe
         procedure      :: destruct           => SVV_destruct
   end type SVV_t

   type(SVV_t), protected :: SVV
   type(Smagorinsky_t)    :: Smagorinsky

   abstract interface
      subroutine Compute_SVV_f(NCONS, NGRAD, Q, Hx, Hy, Hz, sqrt_mu, sqrt_alpha, F)
         use SMConstants, only: RP, NDIM
         implicit none
         integer, intent(in)        :: NCONS, NGRAD
         real(kind=RP), intent(in)  :: Q(NCONS)
         real(kind=RP), intent(in)  :: Hx(NCONS)
         real(kind=RP), intent(in)  :: Hy(NCONS)
         real(kind=RP), intent(in)  :: Hz(NCONS)
         real(kind=RP), intent(in)  :: sqrt_mu
         real(kind=RP), intent(in)  :: sqrt_alpha
         real(kind=RP), intent(out) :: F(NCONS, NDIM)
      end subroutine Compute_SVV_f
      subroutine Compute_Hflux_f(NCONS, NGRAD, Q, Ux, Uy, Uz, sqrt_mu, sqrt_alpha, Hx, Hy, Hz)
         use SMConstants, only: RP, NDIM
         implicit none
         integer,       intent(in)  :: NCONS, NGRAD
         real(kind=RP), intent(in)  :: Q(NCONS), Ux(NGRAD), Uy(NGRAD), Uz(NGRAD)
         real(kind=RP), intent(in)  :: sqrt_mu, sqrt_alpha
         real(kind=RP), intent(out) :: Hx(NCONS), Hy(NCONS), Hz(NCONS)
      end subroutine Compute_Hflux_f
   end interface
!
!  ========
   contains
!  ========
!
      subroutine InitializeSVV(self, controlVariables, mesh)
         use FTValueDictionaryClass
         use mainKeywordsModule
         use PhysicsStorage
         implicit none
         class(SVV_t)                          :: self
         class(FTValueDictionary),  intent(in) :: controlVariables
         class(HexMesh),            intent(in) :: mesh
!
!        ---------------
!        Local variables
!        ---------------
!
         integer     :: eID
         character(len=LINE_LENGTH) :: mu, filter_type
!
!        -------------------------
!        Check if SVV is requested
!        -------------------------
!
         if ( controlVariables % containsKey(SVV_KEY) ) then
            if ( controlVariables % logicalValueForKey(SVV_KEY) ) then
               self % enabled = .true.

            else
               self % enabled = .false.

            end if

         else
            self % enabled = .false.

         end if

         if ( .not. self % enabled ) return
!
!        ---------------------
!        Get the SVV viscosity: the viscosity is later multiplied by 1/N
!        ---------------------
!
         if ( controlVariables % containsKey(SVV_MU1_KEY) ) then

            mu = trim(controlVariables % StringValueForKey(SVV_MU1_KEY,LINE_LENGTH) )
            call ToLower(mu)

            select case ( mu )
               case ('smagorinsky')
                  ! TODO: Use default constructor
                  self % muIsSmagorinsky = .TRUE.
                  self % muSVV = 0.0_RP

                  ! TODO: Use the default constructor
                  Smagorinsky % active = .TRUE.
                  Smagorinsky % requiresWallDistances = .FALSE.
                  Smagorinsky % WallModel = 0  ! No wall model
                  if ( controlVariables % containsKey(SVV_MU_SMAG) ) then
                     Smagorinsky % CS = controlVariables % doublePrecisionValueForKey(SVV_MU_SMAG)
                  else
                     Smagorinsky % CS = 0.2_RP
                  end if

               case default
                  self % muSVV = controlVariables % doublePrecisionValueForKey(SVV_MU_KEY)

            end select

         else
            self % muSVV1 = 0.0_RP

         end if

         if ( controlVariables % containsKey(SVV_MU2_KEY) ) then
            self % muSVV2 = controlVariables % doublePrecisionValueForKey(SVV_MU2_KEY)
         else
            self % muSVV2 = 0.0_RP
         end if

         if ( controlVariables % containsKey(SVV_ALPHA_KEY) ) then
            self % alphaSVV = controlVariables % doublePrecisionValueForKey(SVV_ALPHA_KEY)
         else if ( controlVariables % containsKey(SVV_ALPHA_MU_KEY) ) then
            self % alphaIsPropToMu = .TRUE.
            self % sqrt_mu2alpha   = controlVariables % doublePrecisionValueForKey(SVV_ALPHA_MU_KEY)
            self % alphaSVV        = self % sqrt_mu2alpha * self % muSVV
            self % sqrt_mu2alpha   = sqrt( self % sqrt_mu2alpha )
         else
            self % alphaSVV1 = 0.0_RP
         end if

         if ( controlVariables % containsKey(SVV_ALPHA2_KEY) ) then
            self % alphaSVV2 = controlVariables % doublePrecisionValueForKey(SVV_ALPHA2_KEY)
         else
            self % alphaSVV2 = 0.0_RP
         end if

         self % sqrt_muSVV1    = sqrt(self % muSVV1)
         self % sqrt_muSVV2    = sqrt(self % muSVV2)
         self % sqrt_alphaSVV1 = sqrt(self % alphaSVV1)
         self % sqrt_alphaSVV2 = sqrt(self % alphaSVV2)
!
!        --------------
!        Type of filter
!        --------------
!
         if ( controlVariables % containsKey(FILTER_TYPE_KEY) ) then
            filter_type = controlVariables % StringValueForKey(FILTER_TYPE_KEY,LINE_LENGTH)
            call tolower(filter_type)
            select case ( trim(filter_type) )
               case ("low-pass" ) ; self % filterType = LPASS_FILTER
               case ("high-pass") ; self % filterType = HPASS_FILTER
               case default
                  write(STD_OUT,*) 'ERROR. SVV filter type not recognized. Options are:'
                  write(STD_OUT,*) '   * low-pass'
                  write(STD_OUT,*) '   * high-pass'
                  stop
            end select
         else
            self % filterType = HPASS_FILTER
         end if
!
!        ---------------
!        Shape of filter
!        ---------------
!
         if ( controlVariables % containsKey(FILTER_SHAPE_KEY) ) then
            filter_type = controlVariables % StringValueForKey(FILTER_SHAPE_KEY, LINE_LENGTH)
            call tolower(filter_type)
            select case ( trim(filter_type) )
               case ("power")       ; self % filterShape = POW_FILTER
               case ("sharp")       ; self % filterShape = SHARP_FILTER
               case ("exponential") ; self % filterShape = EXP_FILTER
               case default
                  write(STD_OUT,*) 'ERROR. SVV filter shape not recognized. Options are:'
                  write(STD_OUT,*) '   * power'
                  write(STD_OUT,*) '   * sharp'
                  write(STD_OUT,*) '   * exponential'
                  stop
            end select
         else
            self % filterShape = POW_FILTER
         end if
!
!        -------------------------
!        Get the SVV kernel cutoff
!        -------------------------
!
         if ( controlVariables % containsKey(SVV_CUTOFF_KEY) ) then
            if ( trim(controlVariables % StringValueForKey(SVV_CUTOFF_KEY,LINE_LENGTH)) .eq. "automatic") then
               self % automatic_Psvv = .true.
               self % Psvv = 0.0_RP

            else
               self % automatic_Psvv = .false.
               self % Psvv = controlVariables % doublePrecisionValueForKey(SVV_CUTOFF_KEY)

            end if

         else
            self % automatic_Psvv = .true.
            self % Psvv = 0.0_RP

         end if
!
!        Construct the filters
!        ---------------------
         do eID = 1, mesh % no_of_elements
            associate(Nxyz => mesh % elements(eID) % Nxyz)
            call self % ConstructFilter(Nxyz(1))
            call self % ConstructFilter(Nxyz(2))
            call self % ConstructFilter(Nxyz(3))
            end associate
         end do
!
!        --------------------------------------
!        Select the appropriate HFlux functions
!        --------------------------------------
!
         if ( controlVariables % ContainsKey(DISSIPATION_TYPE_KEY) ) then
            mu = controlVariables % StringValueForKey(DISSIPATION_TYPE_KEY, LINE_LENGTH)
            call tolower(mu)

            select case (trim(mu))
            case ("physical") ; self % diss_type = PHYSICAL_DISS
            case ("guermond") ; self % diss_type = GUERMOND_DISS
            case default
               write(STD_OUT,*) 'ERROR. SVV dissipation type not recognized. Options are:'
               write(STD_OUT,*) '   * Physical'
               write(STD_OUT,*) '   * Guermond'
               errorMessage(STD_OUT)
               stop
            end select
         else
            self % diss_type = PHYSICAL_DISS

         end if

         select case (self % diss_type)
         case (PHYSICAL_DISS)
            select case (grad_vars)
            case(GRADVARS_ENTROPY)
               self % Compute_Hflux => Hflux_physical_dissipation_ENTROPY
               self % Compute_SVV   => SVV_physical_dissipation_ENTROPY
               allocate(self % entropy_indexes(5))
               self % entropy_indexes = [1,2,3,4,5]

            case(GRADVARS_ENERGY)
               self % Compute_Hflux => Hflux_physical_dissipation_ENERGY
               self % Compute_SVV   => SVV_physical_dissipation_ENERGY
               allocate(self % entropy_indexes(3))
               self % entropy_indexes = [2,3,4]

            case default
               write(STD_OUT,*) "ERROR. SVV with physical dissipation is only configured for Energy or Entropy gradient variables"
               errorMessage(STD_OUT)
               stop
            end select

         case (GUERMOND_DISS)
            select case (grad_vars)
            case(GRADVARS_ENTROPY)
               self % Compute_Hflux => Hflux_GuermondPopov_ENTROPY
               self % Compute_SVV   => SVV_GuermondPopov_ENTROPY
               allocate(self % entropy_indexes(5))
               self % entropy_indexes = [1,2,3,4,5]

            case default
               write(STD_OUT,*) "ERROR. SVV with Guermond-Popov (2014) dissipation is only configured for Entropy gradient variables"
               errorMessage(STD_OUT)
               stop
            end select
         end select
!
!        Display the configuration
!        -------------------------
         call self % Describe

      end subroutine InitializeSVV
!
!///////////////////////////////////////////////////////////////////////////////////////////////
!
      subroutine SVV_Describe(this)
         implicit none
         !-arguments----------------------------------------
         class(SVV_t), intent (in) :: this
         !--------------------------------------------------

         if (.not. MPI_Process % isRoot) return

         write(STD_OUT,'(/)')
         call Subsection_Header("Spectral Vanishing Viscosity (SVV)")

         write(STD_OUT,'(30X,A,A30)',advance="no") "->","Dissipation type: "
         select case (this % diss_type)
            case (PHYSICAL_DISS)   ; write(STD_OUT,'(A)') 'Physical'
            case (GUERMOND_DISS)   ; write(STD_OUT,'(A)') 'Guermond'
         end select


         if (this % muIsSmagorinsky) then
            write(STD_OUT,'(30X,A,A30,A,F4.2,A)') "->","Viscosity: ", "Smagorinsky (Cs = ", Smagorinsky % CS,  ")"
            write(STD_OUT,'(30X,A,A30,F10.3)') "->","Alpha viscosity: ", this % alphaSVV1
         else
            write(STD_OUT,'(30X,A,A30,F10.6)') "->","Viscosity: ", this % muSVV
         end if

         if (this % alphaIsPropToMu) then
            write(STD_OUT,'(30X,A,A30,F0.1,A)') "->", "Alpha viscosity: ", (this % sqrt_mu2alpha)**2, "x mu_SVV"
         else
            write(STD_OUT,'(30X,A,A30,F10.6)') "->","Alpha viscosity: ", this % alphaSVV
         end if

         write(STD_OUT,'(30X,A,A30)',advance="no") "->","Filter type: "
         select case (this % filterType)
            case (LPASS_FILTER) ; write(STD_OUT,'(A)') 'low-pass'
            case (HPASS_FILTER) ; write(STD_OUT,'(A)') 'high-pass'
         end select

         write(STD_OUT,'(30X,A,A30)',advance="no") "->","Filter shape: "
         select case (this % filterShape)
            case (POW_FILTER)   ; write(STD_OUT,'(A)') 'power kernel'
            case (EXP_FILTER)   ; write(STD_OUT,'(A)') 'exponential kernel'
            case (SHARP_FILTER) ; write(STD_OUT,'(A)') 'sharp kernel'
         end select

         if (this % automatic_Psvv) then
            write(STD_OUT,'(30X,A,A30,A)') "->","Filter cutoff: ", "automatic"
         else
            write(STD_OUT,'(30X,A,A30,F10.3)') "->","Filter cutoff: ", this % Psvv
         end if

      end subroutine SVV_Describe
!
!///////////////////////////////////////////////////////////////////////////////////////////////
!
      subroutine SVV_ComputeInnerFluxes(self, mesh, e, contravariantFlux)
         use ElementClass
         use PhysicsStorage
         use Physics
         use LESModels
         implicit none
         class(SVV_t)                           :: self
         type(HexMesh)                          :: mesh
         type(Element)                          :: e
         real(kind=RP)           , intent (out) :: contravariantFlux(1:NCONS, 0:e%Nxyz(1), 0:e%Nxyz(2), 0:e%Nxyz(3), 1:NDIM)
!
!        ---------------
!        Local variables
!        ---------------
!
         integer             :: i, j, k, l, fIDs(6), i_f, fID
         real(kind=RP)       :: sqrt_mu(0:e % Nxyz(1), 0:e % Nxyz(2), 0:e % Nxyz(3))
         real(kind=RP)       :: sqrt_alpha(0:e % Nxyz(1), 0:e % Nxyz(2), 0:e % Nxyz(3))
         real(kind=RP)       :: Hx(1:NGRAD, 0:e % Nxyz(1), 0:e % Nxyz(2), 0:e % Nxyz(3))
         real(kind=RP)       :: Hy(1:NGRAD, 0:e % Nxyz(1), 0:e % Nxyz(2), 0:e % Nxyz(3))
         real(kind=RP)       :: Hz(1:NGRAD, 0:e % Nxyz(1), 0:e % Nxyz(2), 0:e % Nxyz(3))
         real(kind=RP)       :: Hxf(1:NGRAD, 0:e % Nxyz(1), 0:e % Nxyz(2), 0:e % Nxyz(3))
         real(kind=RP)       :: Hyf(1:NGRAD, 0:e % Nxyz(1), 0:e % Nxyz(2), 0:e % Nxyz(3))
         real(kind=RP)       :: Hzf(1:NGRAD, 0:e % Nxyz(1), 0:e % Nxyz(2), 0:e % Nxyz(3))
         real(kind=RP)       :: Hxf_aux(1:NGRAD, 0:e % Nxyz(1), 0:e % Nxyz(2), 0:e % Nxyz(3))
         real(kind=RP)       :: Hyf_aux(1:NGRAD, 0:e % Nxyz(1), 0:e % Nxyz(2), 0:e % Nxyz(3))
         real(kind=RP)       :: Hzf_aux(1:NGRAD, 0:e % Nxyz(1), 0:e % Nxyz(2), 0:e % Nxyz(3))
         real(kind=RP)       :: cartesianFlux(1:NCONS, 1:NDIM)
         real(kind=RP)       :: SVV_diss, delta
         real(kind=RP)       :: grad_rho2, rho_sensor, Psvv
         real(kind=RP)       :: Qx(0:e % Nxyz(1),0:e % Nxyz(1))
         real(kind=RP)       :: Qy(0:e % Nxyz(2),0:e % Nxyz(2))
         real(kind=RP)       :: Qz(0:e % Nxyz(3),0:e % Nxyz(3))

!
!        --------------------------------------------
!        Compute the viscosity, using LES if required
!        --------------------------------------------
!
         if ( self % muIsSmagorinsky ) then
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
         else
             sqrt_mu = self % sqrt_muSVV
         end if

         if ( self % alphaIsPropToMu ) then
             sqrt_alpha = self % sqrt_mu2alpha * sqrt_mu
         else
             sqrt_alpha = self % sqrt_alphaSVV
         end if
!
!        ---------------------------------
!        Compute the viscosity and filters
!        ---------------------------------
!
         if ( self % muIsSmagorinsky ) then

            delta = (e % geom % volume / product(e % Nxyz + 1)) ** (1.0_RP/3.0_RP)

            do k = 0, e % Nxyz(3); do j = 0, e % Nxyz(2); do i = 0, e % Nxyz(1)
               call Smagorinsky % ComputeViscosity(delta, e % geom % dWall(i,j,k), &
                                                   e % storage % Q(:,i,j,k),       &
                                                   e % storage % U_x(:,i,j,k),     &
                                                   e % storage % U_y(:,i,j,k),     &
                                                   e % storage % U_z(:,i,j,k),     &
                                                   sqrt_mu(i,j,k))
               sqrt_mu(i,j,k) = sqrt(sqrt_mu(i,j,k))
            end do               ; end do               ; end do

         else
            sqrt_mu = self % sqrt_muSVV1

         end if

         sqrt_alpha = self % sqrt_alphaSVV1

         if ( self % automatic_Psvv ) then
!
!           Use a density sensor to compute the P_svv parameter
!           ---------------------------------------------------
            rho_sensor = 0.0_RP

            do k = 0, e % Nxyz(3); do j = 0, e % Nxyz(2); do i = 0, e % Nxyz(1)

               if (self % diss_type == PHYSICAL_DISS) then
                  grad_rho2 = POW2(e % storage % U_x(1,i,j,k)) &
                            + POW2(e % storage % U_y(1,i,j,k)) &
                            + POW2(e % storage % U_z(1,i,j,k))
               else if (self % diss_type == GUERMOND_DISS) then
                  grad_rho2 = POW2(sum(e % storage % Q(:,i,j,k) * e % storage % U_x(:,i,j,k))) &
                            + POW2(sum(e % storage % Q(:,i,j,k) * e % storage % U_y(:,i,j,k))) &
                            + POW2(sum(e % storage % Q(:,i,j,k) * e % storage % U_z(:,i,j,k)))
               end if

               rho_sensor = rho_sensor + NodalStorage(e % Nxyz(1)) % w(i) &
                                       * NodalStorage(e % Nxyz(2)) % w(j) &
                                       * NodalStorage(e % Nxyz(3)) % w(k) &
                                       * e % geom % jacobian(i,j,k)       &
                                       * grad_rho2

            end do               ; end do               ; end do

            rho_sensor = sqrt(rho_sensor)
!
!           Also update the viscosities if SVV-LES is not active
!           ----------------------------------------------------
            if ( rho_sensor < 1.0_RP ) then
               Psvv = 4.0_RP

               if (.not. self % muIsSmagorinsky ) then
                  sqrt_mu    = self % sqrt_muSVV2
                  sqrt_alpha = self % sqrt_alphaSVV2
               end if

            else
               Psvv = 0.0_RP

               if (.not. self % muIsSmagorinsky ) then
                  sqrt_mu    = self % sqrt_muSVV1
                  sqrt_alpha = self % sqrt_alphaSVV1
               end if

            end if
!
!           Recompute the filters
!           ---------------------
            if ( e % Nxyz(1) > 0 ) then
               call self % filters(e % Nxyz(1)) % Recompute(Psvv, self % filtertype,Qx)
            else
               Qx = self % filters(e % Nxyz(1)) % Q
            end if

            if ( e % Nxyz(2) > 0 ) then
               call self % filters(e % Nxyz(2)) % Recompute(Psvv, self % filtertype,Qy)
            else
               Qy = self % filters(e % Nxyz(2)) % Q
            end if

            if ( e % Nxyz(3) > 0 ) then
               call self % filters(e % Nxyz(3)) % Recompute(Psvv, self % filtertype,Qz)
            else
               Qz = self % filters(e % Nxyz(3)) % Q
            end if

            e % Psvv = Psvv

         else
            e % Psvv = self % Psvv
            Qx = self % filters(e % Nxyz(1)) % Q
            Qy = self % filters(e % Nxyz(2)) % Q
            Qz = self % filters(e % Nxyz(3)) % Q

         end if

         associate(spA_xi   => NodalStorage(e % Nxyz(1)), &
                   spA_eta  => NodalStorage(e % Nxyz(2)), &
                   spA_zeta => NodalStorage(e % Nxyz(3)))
!
!        -----------------
!        Compute the Hflux
!        -----------------
!
         do k = 0, e % Nxyz(3) ; do j = 0, e % Nxyz(2) ; do i = 0, e % Nxyz(1)
           call self % Compute_Hflux(NCONS, NGRAD, e % storage % Q(:,i,j,k), e % storage % U_x(:,i,j,k), &
                                                e % storage % U_y(:,i,j,k), e % storage % U_z(:,i,j,k), &
                                                sqrt_mu(i,j,k), sqrt_alpha(i,j,k), Hx(:,i,j,k), Hy(:,i,j,k), Hz(:,i,j,k))

           Hx(:,i,j,k) = sqrt(e % geom % jacobian(i,j,k)) * Hx(:,i,j,k)
           Hy(:,i,j,k) = sqrt(e % geom % jacobian(i,j,k)) * Hy(:,i,j,k)
           Hz(:,i,j,k) = sqrt(e % geom % jacobian(i,j,k)) * Hz(:,i,j,k)
         end do                ; end do                ; end do
!
!        ----------------
!        Filter the Hflux
!        ----------------
!
!        Perform filtering in xi Hf_aux -> Hf
!        -----------------------
         Hxf_aux = Hx     ; Hyf_aux = Hy     ; Hzf_aux = Hz
         Hxf     = 0.0_RP ; Hyf     = 0.0_RP ; Hzf     = 0.0_RP
         do k = 0, e % Nxyz(3) ; do j = 0, e % Nxyz(2) ; do l = 0, e % Nxyz(1) ; do i = 0, e % Nxyz(1)
               Hxf(:,i,j,k) = Hxf(:,i,j,k) + Qx(i,l) * Hxf_aux(:,l,j,k)
               Hyf(:,i,j,k) = Hyf(:,i,j,k) + Qx(i,l) * Hyf_aux(:,l,j,k)
               Hzf(:,i,j,k) = Hzf(:,i,j,k) + Qx(i,l) * Hzf_aux(:,l,j,k)
         end do                ; end do                ; end do                ; end do
!
!        Perform filtering in eta Hf -> Hf_aux
!        ------------------------
         Hxf_aux = 0.0_RP  ; Hyf_aux = 0.0_RP  ; Hzf_aux = 0.0_RP
         do k = 0, e % Nxyz(3) ; do l = 0, e % Nxyz(2) ; do j = 0, e % Nxyz(2) ; do i = 0, e % Nxyz(1)
            Hxf_aux(:,i,j,k) = Hxf_aux(:,i,j,k) + Qy(j,l) * Hxf(:,i,l,k)
            Hyf_aux(:,i,j,k) = Hyf_aux(:,i,j,k) + Qy(j,l) * Hyf(:,i,l,k)
            Hzf_aux(:,i,j,k) = Hzf_aux(:,i,j,k) + Qy(j,l) * Hzf(:,i,l,k)
         end do                ; end do                ; end do                ; end do
!
!        Perform filtering in zeta Hf_aux -> Hf
!        -------------------------
         Hxf = 0.0_RP  ; Hyf = 0.0_RP  ; Hzf = 0.0_RP
         do l = 0, e % Nxyz(3) ; do k = 0, e % Nxyz(3) ; do j = 0, e % Nxyz(2) ; do i = 0, e % Nxyz(1)
            Hxf(:,i,j,k) = Hxf(:,i,j,k) + Qz(k,l) * Hxf_aux(:,i,j,l)
            Hyf(:,i,j,k) = Hyf(:,i,j,k) + Qz(k,l) * Hyf_aux(:,i,j,l)
            Hzf(:,i,j,k) = Hzf(:,i,j,k) + Qz(k,l) * Hzf_aux(:,i,j,l)
         end do                ; end do                ; end do                ; end do

         if (self % filterType == LPASS_FILTER) then
            Hxf = Hx - Hxf
            Hyf = Hy - Hyf
            Hzf = Hz - Hzf
         end if
!
!        ----------------
!        Compute the flux
!        ----------------
!
         SVV_diss = 0.0_RP
         do k = 0, e%Nxyz(3)   ; do j = 0, e%Nxyz(2) ; do i = 0, e%Nxyz(1)
            call self % Compute_SVV( NCONS, NGRAD, e % storage % Q(:,i,j,k), Hxf(:,i,j,k), Hyf(:,i,j,k), &
                               Hzf(:,i,j,k), sqrt_mu(i,j,k), sqrt_alpha(i,j,k), cartesianFlux )

            cartesianFlux = sqrt(e % geom % invJacobian(i,j,k)) * cartesianFlux
!
!           Scale it with the mesh size
!           ---------------------------
            SVV_diss = SVV_diss + spA_xi % w(i) * spA_eta % w(j) * spA_zeta % w(k) * &
                        (sum(e % storage % U_x(self % entropy_indexes,i,j,k)*cartesianFlux(self % entropy_indexes,IX) + &
                             e % storage % U_y(self % entropy_indexes,i,j,k)*cartesianFlux(self % entropy_indexes,IY) + &
                             e % storage % U_z(self % entropy_indexes,i,j,k)*cartesianFlux(self % entropy_indexes,IZ))) * e % geom % jacobian(i,j,k)



            contravariantFlux(:,i,j,k,IX) =     cartesianFlux(:,IX) * e % geom % jGradXi(IX,i,j,k)  &
                                             +  cartesianFlux(:,IY) * e % geom % jGradXi(IY,i,j,k)  &
                                             +  cartesianFlux(:,IZ) * e % geom % jGradXi(IZ,i,j,k)


            contravariantFlux(:,i,j,k,IY) =     cartesianFlux(:,IX) * e % geom % jGradEta(IX,i,j,k)  &
                                             +  cartesianFlux(:,IY) * e % geom % jGradEta(IY,i,j,k)  &
                                             +  cartesianFlux(:,IZ) * e % geom % jGradEta(IZ,i,j,k)


            contravariantFlux(:,i,j,k,IZ) =     cartesianFlux(:,IX) * e % geom % jGradZeta(IX,i,j,k)  &
                                             +  cartesianFlux(:,IY) * e % geom % jGradZeta(IY,i,j,k)  &
                                             +  cartesianFlux(:,IZ) * e % geom % jGradZeta(IZ,i,j,k)

         end do               ; end do            ; end do

         e % storage % SVV_diss = SVV_diss
!
!        --------------------
!        Prolong to the faces
!        --------------------
!
         fIDs = e % faceIDs
         call e % ProlongHfluxToFaces(NCONS, contravariantFlux, &
                                   mesh % faces(fIDs(1)),&
                                   mesh % faces(fIDs(2)),&
                                   mesh % faces(fIDs(3)),&
                                   mesh % faces(fIDs(4)),&
                                   mesh % faces(fIDs(5)),&
                                   mesh % faces(fIDs(6))   )


         end associate

      end subroutine SVV_ComputeInnerFluxes
!
!//////////////////////////////////////////////////////////////////////////////////////////////////////////////
!
!        Library with Hfluxes
!
!//////////////////////////////////////////////////////////////////////////////////////////////////////////////
!
      subroutine Hflux_physical_dissipation_ENERGY(NCONS, NGRAD, Q, Ux, Uy, Uz, sqrt_mu, sqrt_alpha, Hx, Hy, Hz)
!
!        ***************************************************************************************
!           For the energy variables, the SVV flux is very simple as the NS viscous matrix
!        is constant. We only multiply by the square root of the viscosity
!
!        ***************************************************************************************
!
         implicit none
         integer,    intent(in)     :: NCONS, NGRAD
         real(kind=RP), intent(in)  :: Q(NCONS), Ux(NGRAD), Uy(NGRAD), Uz(NGRAD)
         real(kind=RP), intent(in)  :: sqrt_mu, sqrt_alpha
         real(kind=RP), intent(out) :: Hx(NCONS), Hy(NCONS), Hz(NCONS)
!
!        ---------------
!        Local variables
!        ---------------
!
         real(kind=RP) :: invRho, u, v, w, p_div_rho, sqrt_mu_T

         Hx = sqrt_mu*Ux
         Hy = sqrt_mu*Uy
         Hz = sqrt_mu*Uz

      end subroutine Hflux_physical_dissipation_ENERGY

      subroutine Hflux_physical_dissipation_ENTROPY(NCONS, NGRAD, Q, Ux, Uy, Uz, sqrt_mu, sqrt_alpha, Hx, Hy, Hz)
!
!        ***************************************************************************************
!
!           This Hflux is computed from the LU decomposition of the viscous fluxes.
!        If Fv = Lᵀ·D·L∇U, then Hflux = √D*L∇U, with
!
!        D = diag(α  4/3µT  µT  µT  T²κ | α  0  µT  µT  T²κ | α  0  0  0  T²κ),
!
!        and
!
!            |---------------------|-----------------------|-----------------------|
!            | 1   0   0   0   0   |   0   0   0   0   0   |   0   0   0   0   0   |
!            | 0   1   0   0   u   |   0   0 -1/2  0 -v/2  |   0   0   0 -1/2 -w/2 |
!            | 0   0   1   0   v   |   0   1   0   0   u   |   0   0   0   0   0   |
!            | 0   0   0   1   w   |   0   0   0   0   0   |   0   1   0   0   u   |
!            | 0   0   0   0   1   |   0   0   0   0   0   |   0   0   0   0   0   |
!            |---------------------|-----------------------|-----------------------|
!            | 0   0   0   0   0   |   1   0   0   0   0   |   0   0   0   0   0   |
!            | 0   0   0   0   0   |   0   0   0   0   0   |   0   0   0   0   0   |
!        L = | 0   0   0   0   0   |   0   0   1   0   v   |   0   0   0  -1  -w   |
!            | 0   0   0   0   0   |   0   0   0   1   w   |   0   0   1   0   v   |
!            | 0   0   0   0   0   |   0   0   0   0   1   |   0   0   0   0   0   |
!            |---------------------|-----------------------|-----------------------|
!            | 0   0   0   0   0   |   0   0   0   0   0   |   1   0   0   0   0   |
!            | 0   0   0   0   0   |   0   0   0   0   0   |   0   0   0   0   0   |
!            | 0   0   0   0   0   |   0   0   0   0   0   |   0   0   0   0   0   |
!            | 0   0   0   0   0   |   0   0   0   0   0   |   0   0   0   0   0   |
!            | 0   0   0   0   0   |   0   0   0   0   0   |   0   0   0   0   1   |
!            |---------------------|-----------------------|-----------------------|
!
!        Only the non-constants are taken into the sqrt of (D). (e.g. 4/3µT -> 4/3 √(µT))
!
!        ***************************************************************************************
!
         implicit none
         integer,    intent(in)     :: NCONS, NGRAD
         real(kind=RP), intent(in)  :: Q(NCONS), Ux(NGRAD), Uy(NGRAD), Uz(NGRAD)
         real(kind=RP), intent(in)  :: sqrt_mu, sqrt_alpha
         real(kind=RP), intent(out) :: Hx(NCONS), Hy(NCONS), Hz(NCONS)
!
!        ---------------
!        Local variables
!        ---------------
!
         real(kind=RP) :: invRho, u, v, w, p_div_rho, sqrt_mu_T, mu_to_kappa_gammaM2

         invRho = 1.0_RP / Q(IRHO)

         u = Q(IRHOU) * invRho
         v = Q(IRHOV) * invRho
         w = Q(IRHOW) * invRho

         p_div_rho = thermodynamics % gammaMinus1*(invRho * Q(IRHOE)-0.5_RP*(u*u+v*v+w*w))
         sqrt_mu_T = sqrt_mu*sqrt(p_div_rho)

         mu_to_kappa_gammaM2 = dimensionless % mu_to_kappa * dimensionless % gammaM2

         Hx(IRHO)  = Ux(IRHO)
         Hx(IRHOU) = Ux(IRHOU) + u*Ux(IRHOE) - 0.5_RP*(Uy(IRHOV) + v*Uy(IRHOE) + Uz(IRHOW) + w*Uz(IRHOE))
         Hx(IRHOV) = Ux(IRHOV) + v*Ux(IRHOE) + Uy(IRHOU) + u*Uy(IRHOE)
         Hx(IRHOW) = Ux(IRHOW) + w*Ux(IRHOE) + Uz(IRHOU) + u*Uz(IRHOE)
         Hx(IRHOE) = Ux(IRHOE)

         Hy(IRHO)  = Uy(IRHO)
         Hy(IRHOU) = 0.0_RP
         Hy(IRHOV) = Uy(IRHOV) + v*Uy(IRHOE) - Uz(IRHOW) - w*Uz(IRHOE)
         Hy(IRHOW) = Uy(IRHOW) + w*Uy(IRHOE) + Uz(IRHOV) + v*Uz(IRHOE)
         Hy(IRHOE) = Uy(IRHOE)

         Hz(IRHO)  = Uz(IRHO)
         Hz(IRHOU) = 0.0_RP
         Hz(IRHOV) = 0.0_RP
         Hz(IRHOW) = 0.0_RP
         Hz(IRHOE) = Uz(IRHOE)

         Hx = Hx*[sqrt_alpha, 4.0_RP/3.0_RP*sqrt_mu_T, sqrt_mu_T, sqrt_mu_T, sqrt_mu*mu_to_kappa_gammaM2*p_div_rho]
         Hy = Hy*[sqrt_alpha, 0.0_RP                 , sqrt_mu_T, sqrt_mu_T, sqrt_mu*mu_to_kappa_gammaM2*p_div_rho]

         Hz(IRHO)  = Hz(IRHO)*sqrt_alpha
         Hz(IRHOE) = Hz(IRHOE)*sqrt_mu*mu_to_kappa_gammaM2*p_div_rho

      end subroutine Hflux_physical_dissipation_ENTROPY

      subroutine Hflux_GuermondPopov_ENTROPY(NCONS, NGRAD, Q, Ux, Uy, Uz, sqrt_mu, sqrt_alpha, Hx, Hy, Hz)
!
!        ***************************************************************************************
!
!           This Hflux is computed from the LU decomposition of the Guermond-Popov fluxes.
!        If FGP = Lᵀ·D·L∇U, then Hflux = √D*L∇U, with
!
!        D = diag(αρ  µp  µp/2  µp/2  αρ | αρ  0  µp  µp/2  αρ  | αρ  0  0  µp  αρ),
!
!        and
!
!            |---------------------|-----------------------|-----------------------|
!            | 1   u   v   w   e   |   0   0   0   0   0   |   0   0   0   0   0   |
!            | 0   1   0   0   u   |   0   0   0   0   0   |   0   0   0   0   0   |
!            | 0   0   1   0   v   |   0   1   0   0   u   |   0   0   0   0   0   |
!            | 0   0   0   1   w   |   0   0   0   0   0   |   0   1   0   0   u   |
!            | 0   0   0   0   Λ   |   0   0   0   0   0   |   0   0   0   0   0   |
!            |---------------------|-----------------------|-----------------------|
!            | 0   0   0   0   0   |   1   u   v   w   e   |   0   0   0   0   0   |
!            | 0   0   0   0   0   |   0   0   0   0   0   |   0   0   0   0   0   |
!        L = | 0   0   0   0   0   |   0   0   1   0   v   |   0   0   0   0   0   |
!            | 0   0   0   0   0   |   0   0   0   1   w   |   0   0   1   0   v   |
!            | 0   0   0   0   0   |   0   0   0   0   Λ   |   0   0   0   0   0   |
!            |---------------------|-----------------------|-----------------------|
!            | 0   0   0   0   0   |   0   0   0   0   0   |   1   u   v   w   e   |
!            | 0   0   0   0   0   |   0   0   0   0   0   |   0   0   0   0   0   |
!            | 0   0   0   0   0   |   0   0   0   0   0   |   0   0   0   0   0   |
!            | 0   0   0   0   0   |   0   0   0   0   0   |   0   0   0   1   w   |
!            | 0   0   0   0   0   |   0   0   0   0   0   |   0   0   0   0   Λ   |
!            |---------------------|-----------------------|-----------------------|
!
!        Λ=(p/ρ)/√(γ-1)
!
!        Only the non-constants are taken into the sqrt of (D). (e.g. µp/2 -> 1/2 √(µp) )
!
!        ***************************************************************************************
!
         implicit none
         integer,    intent(in)     :: NCONS, NGRAD
         real(kind=RP), intent(in)  :: Q(NCONS), Ux(NGRAD), Uy(NGRAD), Uz(NGRAD)
         real(kind=RP), intent(in)  :: sqrt_mu, sqrt_alpha
         real(kind=RP), intent(out) :: Hx(NCONS), Hy(NCONS), Hz(NCONS)
!
!        ---------------
!        Local variables
!        ---------------
!
         real(kind=RP) :: invRho, u, v, w, p_div_rho, lambda
         real(kind=RP), parameter :: inv_sqrt_gamma_minus_1 = 1.0_RP / sqrt(1.4_RP-1.0_RP)
         real(kind=RP) :: sqrt_alpha_rho, sqrt_mu_p

         invRho = 1.0_RP / Q(IRHO)

         u = Q(IRHOU) * invRho
         v = Q(IRHOV) * invRho
         w = Q(IRHOW) * invRho

         p_div_rho = thermodynamics % gammaMinus1*(invRho * Q(IRHOE)-0.5_RP*(u*u+v*v+w*w))
         lambda = inv_sqrt_gamma_minus_1*p_div_rho

         sqrt_alpha_rho = sqrt_alpha*sqrt(Q(IRHO))
         sqrt_mu_p      = sqrt_mu*sqrt(p_div_rho*Q(IRHO))

         Hx(IRHO)  = invRho*sum(Q*Ux)
         Hx(IRHOU) = Ux(IRHOU) + u*Ux(IRHOE)
         Hx(IRHOV) = Ux(IRHOV) + v*Ux(IRHOE) + Uy(IRHOU) + u*Uy(IRHOE)
         Hx(IRHOW) = Ux(IRHOW) + w*Ux(IRHOE) + Uz(IRHOU) + u*Uz(IRHOE)
         Hx(IRHOE) = lambda*Ux(IRHOE)

         Hy(IRHO)  = invRho*sum(Q*Uy)
         Hy(IRHOU) = 0.0_RP
         Hy(IRHOV) = Uy(IRHOV) + v*Uy(IRHOE)
         Hy(IRHOW) = Uy(IRHOW) + w*Uy(IRHOE) + Uz(IRHOV) + v*Uz(IRHOE)
         Hy(IRHOE) = lambda*Uy(IRHOE)

         Hz(IRHO)  = invRho*sum(Q*Uz)
         Hz(IRHOU) = 0.0_RP
         Hz(IRHOV) = 0.0_RP
         Hz(IRHOW) = Uz(IRHOW) + w*Uz(IRHOE)
         Hz(IRHOE) = lambda*Uz(IRHOE)
!
!        Scale with sqrt(D)
!        ------------------
         Hx = Hx*[sqrt_alpha_rho, sqrt_mu_p, 0.5_RP*sqrt_mu_p, 0.5_RP*sqrt_mu_p, sqrt_alpha_rho]
         Hy = Hy*[sqrt_alpha_rho, 0.0_RP   , sqrt_mu_p       , 0.5_RP*sqrt_mu_p, sqrt_alpha_rho]
         Hz = Hz*[sqrt_alpha_rho, 0.0_RP   , 0.0_RP          , sqrt_mu_p       , sqrt_alpha_rho]

      end subroutine Hflux_GuermondPopov_ENTROPY

!
!//////////////////////////////////////////////////////////////////////////////////////////////////////////////
!
!        Library with SVV dissipations f(Q,H)
!
!//////////////////////////////////////////////////////////////////////////////////////////////////////////////
!
      subroutine SVV_physical_dissipation_ENERGY(NCONS, NGRAD, Q, Hx, Hy, Hz, sqrt_mu, sqrt_alpha, F)
         implicit none
         integer, intent(in)        :: NCONS, NGRAD
         real(kind=RP), intent(in)  :: Q(NCONS)
         real(kind=RP), intent(in)  :: Hx(NCONS)
         real(kind=RP), intent(in)  :: Hy(NCONS)
         real(kind=RP), intent(in)  :: Hz(NCONS)
         real(kind=RP), intent(in)  :: sqrt_mu
         real(kind=RP), intent(in)  :: sqrt_alpha
         real(kind=RP), intent(out) :: F(NCONS, NDIM)
!
!        ---------------
!        Local variables
!        ---------------
!
         real(kind=RP)  :: invRho, divV, u(2:1+NDIM)
         real(kind=RP)  :: kappa

         kappa = sqrt_mu * dimensionless % mu_to_kappa

         invRho  = 1.0_RP / Q(IRHO)
         u = Q(IRHOU:IRHOW)*invRho

         divV = Hx(IX) + Hy(IY) + Hz(IZ)

         F(IRHO,IX)  = 0.0_RP
         F(IRHOU,IX) = sqrt_mu * (2.0_RP * Hx(IRHOU) - 2.0_RP/3.0_RP * divV )
         F(IRHOV,IX) = sqrt_mu * ( Hx(IRHOV) + Hy(IRHOU) )
         F(IRHOW,IX) = sqrt_mu * ( Hx(IRHOW) + Hz(IRHOU) )
         F(IRHOE,IX) = F(IRHOU,IX) * u(1) + F(IRHOV,IX) * u(2) + F(IRHOW,IX) * u(3) + kappa * Hx(IRHOE)

         F(IRHO,IY) = 0.0_RP
         F(IRHOU,IY) = F(IRHOV,IX)
         F(IRHOV,IY) = sqrt_mu * (2.0_RP * Hy(IRHOV) - 2.0_RP / 3.0_RP * divV )
         F(IRHOW,IY) = sqrt_mu * ( Hy(IRHOW) + Hz(IRHOV) )
         F(IRHOE,IY) = F(IRHOU,IY) * u(1) + F(IRHOV,IY) * u(2) + F(IRHOW,IY) * u(3) + kappa * Hy(IRHOE)

         F(IRHO,IZ) = 0.0_RP
         F(IRHOU,IZ) = F(IRHOW,IX)
         F(IRHOV,IZ) = F(IRHOW,IY)
         F(IRHOW,IZ) = sqrt_mu * ( 2.0_RP * Hz(IRHOW) - 2.0_RP / 3.0_RP * divV )
         F(IRHOE,IZ) = F(IRHOU,IZ) * u(1) + F(IRHOV,IZ) * u(2) + F(IRHOW,IZ) * u(3) + kappa * Hz(IRHOE)

      end subroutine SVV_physical_dissipation_ENERGY

      subroutine SVV_physical_dissipation_ENTROPY(NCONS, NGRAD, Q, Hx, Hy, Hz, sqrt_mu, sqrt_alpha, F)
!
!        ***************************************************************************************
!
!           We add what remains from the decomposition, from Hflux: Fv = Lᵀ·√D·H.
!
!        Recall that in D we took away the constants (to avoid innecesary sqrts)
!
!        D = diag(α  µT  µT  µT  T²κ | α  0  µT  µT  T²κ | α  0  0  0  T²κ),
!
!        and
!
!            |---------------------|-----------------------|-----------------------|
!            | 1   0   0   0   0   |   0   0   0   0   0   |   0   0   0   0   0   |
!            | 0   1   0   0   0   |   0   0   0   0   0   |   0   0   0   0   0   |
!            | 0   0   1   0   0   |   0   0   0   0   0   |   0   0   0   0   0   |
!            | 0   0   0   1   0   |   0   0   0   0   0   |   0   0   0   0   0   |
!            | 0   u   v   w   1   |   0   0   0   0   0   |   0   0   0   0   0   |
!            |---------------------|-----------------------|-----------------------|
!            | 0   0   0   0   0   |   1   0   0   0   0   |   0   0   0   0   0   |
!            | 0   0   1   0   0   |   0   0   0   0   0   |   0   0   0   0   0   |
!        Lᵀ= | 0 -1/2  0   0   0   |   0   0   1   0   0   |   0   0   0   0   0   |
!            | 0   0   0   0   0   |   0   0   0   1   0   |   0   0   0   0   0   |
!            | 0 -v/2  u   0   0   |   0   0   v   w   1   |   0   0   0   0   0   |
!            |---------------------|-----------------------|-----------------------|
!            | 0   0   0   0   0   |   0   0   0   0   0   |   1   0   0   0   0   |
!            | 0   0   0   1   0   |   0   0   0   0   0   |   0   0   0   0   0   |
!            | 0   0   0   0   0   |   0   0   0   1   0   |   0   0   0   0   0   |
!            | 0 -1/2  0   0   0   |   0   0  -1   0   0   |   0   0   0   0   0   |
!            | 0 -w/2  0   u   0   |   0   0  -w   v   0   |   0   0   0   0   1   |
!            |---------------------|-----------------------|-----------------------|
!
!        ***************************************************************************************
!

         implicit none
         integer, intent(in)        :: NCONS, NGRAD
         real(kind=RP), intent(in)  :: Q(NCONS)
         real(kind=RP), intent(in)  :: Hx(NCONS)
         real(kind=RP), intent(in)  :: Hy(NCONS)
         real(kind=RP), intent(in)  :: Hz(NCONS)
         real(kind=RP), intent(in)  :: sqrt_mu
         real(kind=RP), intent(in)  :: sqrt_alpha
         real(kind=RP), intent(out) :: F(NCONS, NDIM)
!
!        ---------------
!        Local variables
!        ---------------
!
         real(kind=RP) :: Hx_sqrtD(NCONS), Hy_sqrtD(NCONS), Hz_sqrtD(NCONS)
         real(kind=RP) :: invRho, u, v, w, p_div_rho, sqrt_mu_T

         invRho = 1.0_RP / Q(IRHO)

         u = Q(IRHOU) * invRho
         v = Q(IRHOV) * invRho
         w = Q(IRHOW) * invRho

         p_div_rho = thermodynamics % gammaMinus1*(invRho * Q(IRHOE)-0.5_RP*(u*u+v*v+w*w))
         sqrt_mu_T = sqrt_mu*sqrt(p_div_rho)

         Hx_sqrtD = Hx*[sqrt_alpha, sqrt_mu_T, sqrt_mu_T, sqrt_mu_T, sqrt_mu*p_div_rho]
         Hy_sqrtD = Hy*[sqrt_alpha, 0.0_RP   , sqrt_mu_T, sqrt_mu_T, sqrt_mu*p_div_rho]

         Hz_sqrtD(IRHO)        = Hz(IRHO)*sqrt_alpha
         Hz_sqrtD(IRHOU:IRHOW) = 0.0_RP
         Hz_sqrtD(IRHOE)       = Hz(IRHOE)*sqrt_mu*p_div_rho

         F(IRHO,IX)  = Hx_sqrtD(IRHO)
         F(IRHOU,IX) = Hx_sqrtD(IRHOU)
         F(IRHOV,IX) = Hx_sqrtD(IRHOV)
         F(IRHOW,IX) = Hx_sqrtD(IRHOW)
         F(IRHOE,IX) = u*Hx_sqrtD(IRHOU) + v*Hx_sqrtD(IRHOV) + w*Hx_sqrtD(IRHOW) + Hx_sqrtD(IRHOE)

         F(IRHO,IY)  = Hy_sqrtD(IRHO)
         F(IRHOU,IY) = Hx_sqrtD(IRHOV)
         F(IRHOV,IY) = -0.5_RP*Hx_sqrtD(IRHOU)+Hy_sqrtD(IRHOV)
         F(IRHOW,IY) = Hy_sqrtD(IRHOW)
         F(IRHOE,IY) = -0.5_RP*v*Hx_sqrtD(IRHOU) + u*Hx_sqrtD(IRHOV) + v*Hy_sqrtD(IRHOV) + w*Hy_sqrtD(IRHOW) + Hy_sqrtD(IRHOE)

         F(IRHO,IZ)  = Hz_sqrtD(IRHO)
         F(IRHOU,IZ) = Hx_sqrtD(IRHOW)
         F(IRHOV,IZ) = Hy_sqrtD(IRHOW)
         F(IRHOW,IZ) = -0.5_RP*Hx_sqrtD(IRHOU) - Hy_sqrtD(IRHOV)
         F(IRHOE,IZ) = -0.5_RP*w*Hx_sqrtD(IRHOU) + u*Hx_sqrtD(IRHOW) - w*Hy_sqrtD(IRHOV) + v*Hy_sqrtD(IRHOW) + Hz_sqrtD(IRHOE)

      end subroutine SVV_physical_dissipation_ENTROPY

      subroutine SVV_GuermondPopov_ENTROPY(NCONS, NGRAD, Q, Hx, Hy, Hz, sqrt_mu, sqrt_alpha, F)
!
!        ***************************************************************************************
!
!           We add what remains from the decomposition, from Hflux: FGP = Lᵀ·√D·H.
!
!        Recall that in D we took away the constants (to avoid innecesary sqrts)
!
!        D = diag(αρ  µp  µp  µp  αρ | αρ  0  µp  µp  αρ  | αρ  0  0  µp  αρ),
!
!        and
!
!            |---------------------|-----------------------|-----------------------|
!            | 1   0   0   0   0   |   0   0   0   0   0   |   0   0   0   0   0   |
!            | u   1   0   0   0   |   0   0   0   0   0   |   0   0   0   0   0   |
!            | v   0   1   0   0   |   0   0   0   0   0   |   0   0   0   0   0   |
!            | w   0   0   1   0   |   0   0   0   0   0   |   0   0   0   0   0   |
!            | e   u   v   w   Λ   |   0   0   0   0   0   |   0   0   0   0   0   |
!            |---------------------|-----------------------|-----------------------|
!            | 0   0   0   0   0   |   1   0   0   0   0   |   0   0   0   0   0   |
!            | 0   0   1   0   0   |   u   0   0   0   0   |   0   0   0   0   0   |
!        Lᵀ= | 0   0   0   0   0   |   v   0   1   0   0   |   0   0   0   0   0   |
!            | 0   0   0   0   0   |   w   0   0   1   0   |   0   0   0   0   0   |
!            | 0   0   u   0   0   |   e   0   v   w   Λ   |   0   0   0   0   0   |
!            |---------------------|-----------------------|-----------------------|
!            | 0   0   0   0   0   |   0   0   0   0   0   |   1   0   0   0   0   |
!            | 0   0   0   1   0   |   0   0   0   0   0   |   u   0   0   0   0   |
!            | 0   0   0   0   0   |   0   0   0   1   0   |   v   0   0   0   0   |
!            | 0   0   0   0   0   |   0   0   0   0   0   |   w   0   0   1   0   |
!            | 0   0   0   u   0   |   0   0   0   v   0   |   e   0   0   w   Λ   |
!            |---------------------|-----------------------|-----------------------|
!
!        ***************************************************************************************
!

         implicit none
         integer, intent(in)        :: NCONS, NGRAD
         real(kind=RP), intent(in)  :: Q(NCONS)
         real(kind=RP), intent(in)  :: Hx(NCONS)
         real(kind=RP), intent(in)  :: Hy(NCONS)
         real(kind=RP), intent(in)  :: Hz(NCONS)
         real(kind=RP), intent(in)  :: sqrt_mu
         real(kind=RP), intent(in)  :: sqrt_alpha
         real(kind=RP), intent(out) :: F(NCONS, NDIM)
!
!        ---------------
!        Local variables
!        ---------------
!
         real(kind=RP)            :: Hx_sqrtD(NCONS), Hy_sqrtD(NCONS), Hz_sqrtD(NCONS)
         real(kind=RP)            :: invRho, u, v, w, p_div_rho, lambda, e
         real(kind=RP)            :: sqrt_alpha_rho, sqrt_mu_p
         real(kind=RP), parameter :: inv_sqrt_gamma_minus_1 = 1.0_RP / sqrt(1.4_RP-1.0_RP)

         invRho = 1.0_RP / Q(IRHO)

         u = Q(IRHOU) * invRho
         v = Q(IRHOV) * invRho
         w = Q(IRHOW) * invRho
         e = Q(IRHOE) * invRho

         p_div_rho = thermodynamics % gammaMinus1*(invRho * Q(IRHOE)-0.5_RP*(u*u+v*v+w*w))
         lambda = inv_sqrt_gamma_minus_1*p_div_rho

         sqrt_alpha_rho = sqrt_alpha*sqrt(Q(IRHO))
         sqrt_mu_p      = sqrt_mu*sqrt(p_div_rho*Q(IRHO))

         Hx_sqrtD = Hx*[sqrt_alpha_rho, sqrt_mu_p, sqrt_mu_p, sqrt_mu_p, sqrt_alpha_rho]
         Hy_sqrtD = Hy*[sqrt_alpha_rho, 0.0_RP   , sqrt_mu_p, sqrt_mu_p, sqrt_alpha_rho]
         Hz_sqrtD = Hz*[sqrt_alpha_rho, 0.0_RP   , 0.0_RP   , sqrt_mu_p, sqrt_alpha_rho]

         F(IRHO,IX)  = Hx_sqrtD(IRHO)
         F(IRHOU,IX) = u*Hx_sqrtD(IRHO) + Hx_sqrtD(IRHOU)
         F(IRHOV,IX) = v*Hx_sqrtD(IRHO) + Hx_sqrtD(IRHOV)
         F(IRHOW,IX) = w*Hx_sqrtD(IRHO) + Hx_sqrtD(IRHOW)
         F(IRHOE,IX) = e*Hx_sqrtD(IRHO) + u*Hx_sqrtD(IRHOU) + v*Hx_sqrtD(IRHOV) + w*Hx_sqrtD(IRHOW) + lambda*Hx_sqrtD(IRHOE)

         F(IRHO,IY)  = Hy_sqrtD(IRHO)
         F(IRHOU,IY) = u*Hy_sqrtD(IRHO) + Hx_sqrtD(IRHOV)
         F(IRHOV,IY) = v*Hy_sqrtD(IRHO) + Hy_sqrtD(IRHOV)
         F(IRHOW,IY) = w*Hy_sqrtD(IRHO) + Hy_sqrtD(IRHOW)
         F(IRHOE,IY) = e*Hy_sqrtD(IRHO) + u*Hx_sqrtD(IRHOV) + v*Hy_sqrtD(IRHOV) + w*Hy_sqrtD(IRHOW) + lambda*Hy_sqrtD(IRHOE)

         F(IRHO,IZ)  = Hz_sqrtD(IRHO)
         F(IRHOU,IZ) = u*Hz_sqrtD(IRHO) + Hx_sqrtD(IRHOW)
         F(IRHOV,IZ) = v*Hz_sqrtD(IRHO) + Hy_sqrtD(IRHOW)
         F(IRHOW,IZ) = w*Hz_sqrtD(IRHO) + Hz_sqrtD(IRHOW)
         F(IRHOE,IZ) = e*Hz_sqrtD(IRHO) + u*Hx_sqrtD(IRHOW) + v*Hy_sqrtD(IRHOW) + w*Hz_sqrtD(IRHOW) + lambda*Hz_sqrtD(IRHOE)

      end subroutine SVV_GuermondPopov_ENTROPY
!
!//////////////////////////////////////////////////////////////////////////////////////////////////
!
!           Auxiliar subroutines
!           --------------------
!
!//////////////////////////////////////////////////////////////////////////////////////////////////
!
      subroutine SVV_constructFilter(self, N)
         implicit none
         class(SVV_t)         :: self
         integer, intent(in) :: N
!
!        ---------------
!        Local variables
!        ---------------
!
         integer        :: i, j, k
         integer        :: sharpCutOff
         real(kind=RP)  :: Nodal2Modal(0:N,0:N)
         real(kind=RP)  :: Modal2Nodal(0:N,0:N)
         real(kind=RP)  :: filterCoefficients(0:N)
         real(kind=RP)  :: Lkj(0:N,0:N), dLk_dummy
         real(kind=RP)  :: normLk(0:N)

         if ( self % filters(N) % Constructed ) return

         if ( N .eq. 0 ) then
            self % filters(N) % N = N
            allocate(self % filters(N) % Q(0:N,0:N))
            allocate(self % filters(N) % F(0:N,0:N))
            allocate(self % filters(N) % B(0:N,0:N))
            self % filters(N) % Q = 1.0_RP
            self % filters(N) % F = 1.0_RP
            self % filters(N) % B = 1.0_RP
            self % filters(N) % constructed = .true.
            return
         end if

!
!        Get the evaluation of Legendre polynomials at the interpolation nodes
!        ---------------------------------------------------------------------
         do j = 0 , N ;    do k = 0 , N
            call LegendrePolyAndDerivative(k, NodalStorage(N) % x(j), Lkj(k,j), dLk_dummy)
         end do       ;    end do
!
!        Get the norm of Legendre polynomials
!        ------------------------------------
         normLk = 0.0_RP
         do k = 0 , N   ; do j = 0 , N
            normLk(k) = normLk(k) + NodalStorage(N) % w(j) * Lkj(k,j) * Lkj(k,j)
         end do         ; end do
!
!        Get the transformation from Nodal to Modal and viceversa matrices
!        -----------------------------------------------------------------
         do k = 0 , N   ; do i = 0 , N
            Nodal2Modal(k,i) = NodalStorage(N) % w(i) * Lkj(k,i) / normLk(k)
            Modal2Nodal(i,k) = Lkj(k,i)
         end do         ; end do
!
!        Get the filter coefficients
!        ---------------------------
         select case (self % filterShape)
            case (POW_FILTER)
               do k = 0, N
                  filterCoefficients(k) = (real(k, kind=RP) / N + 1.0e-12_RP) ** self % Psvv
               end do

            case (SHARP_FILTER)
               sharpCutOff = nint(self % Psvv)

               filterCoefficients = 0._RP

               do k = 0, N
                  if ( k >= sharpCutOff ) filterCoefficients(k) = 1.0_RP
               end do

            case (EXP_FILTER)
               filterCoefficients = 0._RP
               do k = 0, N
                  if (k > self % Psvv) filterCoefficients(k) = exp( -real( (k-N)**2 , kind=RP) / (k - self % Psvv) ** 2 )
               end do

         end select

         if (self % filterType == LPASS_FILTER) then
            filterCoefficients = 1._RP - filterCoefficients
         end if
!
!        Compute the filtering matrix
!        ----------------------------
         self % filters(N) % N = N
         allocate(self % filters(N) % Q(0:N,0:N))
         allocate(self % filters(N) % F(0:N,0:N))
         allocate(self % filters(N) % B(0:N,0:N))

         self % filters(N) % Q = 0.0_RP
         do k = 0, N ; do j = 0, N  ; do i = 0, N
            self % filters(N) % Q(i,j) = self % filters(N) % Q(i,j) + Modal2Nodal(i,k) * filterCoefficients(k) * Nodal2Modal(k,j)
         end do      ; end do       ; end do

         self % filters(N) % F = Nodal2Modal
         self % filters(N) % B = Modal2Nodal

         self % filters(N) % constructed = .true.

      end subroutine SVV_constructFilter

      subroutine SVV_destruct(this)
         implicit none
         class(SVV_t) :: this
         integer :: i

         do i = 0, Nmax
            if ( this % filters(i) % constructed ) deallocate(this % filters(i) % Q)
         end do

         !if (this % muIsSmagorinsky) call Smagorinsky % destruct

      end subroutine SVV_destruct

      subroutine FilterMatrices_Recompute(self,Psvv, type_, Q)
         implicit none
         class(FilterMatrices_t)   :: self
         real(kind=RP), intent(in) :: Psvv
         integer,       intent(in) :: type_
         real(kind=RP), intent(out) :: Q(0:self % N,0:self % N)
         integer :: i,j,k
         real(kind=RP)  :: filterCoefficients(0:self % N)

         if (Psvv > 100.0_RP ) then
            filterCoefficients = 0.0_RP
         else
            do k = 0, self % N
               filterCoefficients(k) = (real(k, kind=RP) / self % N + 1.0e-12_RP) ** Psvv
            end do
         end if

         if ( type_ == LPASS_FILTER) then
            filterCoefficients = 1._RP - filterCoefficients
         end if

         Q = 0.0_RP
         do k = 0, self % N ; do j = 0, self % N  ; do i = 0, self % N
            Q(i,j) = Q(i,j) + self % B(i,k) * filterCoefficients(k) * self % F(k,j)
         end do      ; end do       ; end do

      end subroutine FilterMatrices_Recompute

end module SpectralVanishingViscosity
#endif
