#include "Includes.h"
module LESModels
   use SMConstants
   use PhysicsStorage_NS
   use FTValueDictionaryClass
   use Physics_NSKeywordsModule
   use MPI_Process_Info
   use Headers
   use Utilities                 , only: toLower
   use FluidData_NS
   use VariableConversion_NS     , only: getVelocityGradients, getTemperatureGradient, getVelocityGradients_State
   implicit none

   private
   public LESModel, InitializeLESModel, Smagorinsky_t, LESModel_Selector, Smagorinsky_ComputeViscosity
   
!  Keywords
!  --------
   character(len=*), parameter   :: LESIntensityKey = "les model intensity"
   character(len=*), parameter   :: WallModelKey    = "wall model"
   
!  Model parameters
!  ----------------
   real(kind=RP)   , parameter   :: K_VONKARMAN     = 0.4_RP
   !$acc declare copyin(K_VONKARMAN)
   
!  Wall models
!  -----------
   integer         , parameter   :: NO_WALLMODEL       = 0
   integer         , parameter   :: LINEAR_WALLMODEL   = 1
   
   type LESModel_t
      logical  :: active
      logical  :: requiresWallDistances = .FALSE.
      integer  :: WallModel
      integer  :: LESModel
      real(kind=RP)  :: C
      contains
         procedure            :: Initialize         => LESModel_Initialize
         ! procedure, private   :: ComputeWallEffect  => LESModel_ComputeWallEffect
         ! procedure            :: ComputeViscosity   => LESModel_ComputeViscosity
         procedure            :: Describe           => LESModel_Describe
   end type LESModel_t

   type, extends(LESModel_t)  :: Smagorinsky_t
      ! real(kind=RP)  :: CS
      contains
         procedure          :: Initialize         => Smagorinsky_Initialize
         procedure          :: Describe           => Smagorinsky_Describe
         ! procedure          :: ComputeViscosity   => Smagorinsky_ComputeViscosity
   end type Smagorinsky_t

   type, extends(LESModel_t)  :: WALE_t
      ! real(kind=RP)  :: Cw
      contains
         procedure          :: Initialize         => WALE_Initialize
         procedure          :: Describe           => WALE_Describe
         ! procedure          :: ComputeViscosity   => WALE_ComputeViscosity
   end type WALE_t

   type, extends(LESModel_t)  :: Vreman_t
      ! real(kind=RP)  :: C
      contains
         procedure          :: Initialize         => Vreman_Initialize
         procedure          :: Describe           => Vreman_Describe
         ! procedure          :: ComputeViscosity   => Vreman_ComputeViscosity
   end type Vreman_t


   class(LESModel_t), allocatable   :: LESModel

   contains
      subroutine InitializeLESModel(model, controlVariables)
         implicit none
         class(LESModel_t), allocatable        :: model
         class(FTValueDictionary),  intent(in) :: controlVariables
!
!        ---------------
!        Local variables         
!        ---------------
!
         character(len=LINE_LENGTH)    :: modelName

!
!        Select LES model
!        ----------------
         if ( controlVariables % containsKey(LESMODEL_KEY) ) then
            modelName = controlVariables % stringValueForKey(LESMODEL_KEY, LINE_LENGTH)
            call toLower(modelName)

            select case (trim(modelName))
            case ("none")
               if (.not. allocated(model)) allocate(LESModel_t     :: model)
               model % LESModel = 0 
            case ("smagorinsky")
               if (.not. allocated(model)) allocate(Smagorinsky_t  :: model)
               model % LESModel = 1 

            case ("wale")
               if (.not. allocated(model)) allocate(WALE_t  :: model)
               model % LESModel = 2 

            case ("vreman")
               if (.not. allocated(model)) allocate(Vreman_t  :: model)
               model % LESModel = 3 

            case default
               write(STD_OUT,'(A,A,A)') "LES Model ",trim(modelName), " is not implemented."
               print*, "Available options are:"
               print*, "   * None (default)"
               print*, "   * Smagorinsky"
               print*, "   * Wale"
               print*, "   * Vreman"
               errorMessage(STD_OUT)
               error stop

            end select
   
         else
            if (.not. allocated(model)) allocate(LESModel_t  :: model)

         end if

         call model % Initialize(controlVariables)
         
!        Select wall model
!        -----------------
         if ( controlVariables % containsKey(WallModelKey) ) then
            modelName = controlVariables % stringValueForKey(WallModelKey, LINE_LENGTH)
            call toLower(modelName)

            select case (trim(modelName))
            case ("none")
               model % WallModel = NO_WALLMODEL
               model % requiresWallDistances = .true.

            case ("linear")
               model % WallModel             = LINEAR_WALLMODEL
               model % requiresWallDistances = .true.

            case ("Wale")
               model % WallModel = NO_WALLMODEL
               model % requiresWallDistances = .true.

            case ("Vreman")
               model % WallModel = NO_WALLMODEL
               model % requiresWallDistances = .true.

            case default
               write(STD_OUT,'(A,A,A)') "Wall model ",trim(modelName), " is not implemented."
               print*, "Available options are:"
               print*, "   * Linear"
               print*, "   * Wale"
               print*, "   * Vreman"
               print*, "   * None (default)"
               errorMessage(STD_OUT)
               error stop

            end select
            
         else
            !model % WallModel = LINEAR_WALLMODEL
            model % WallModel = NO_WALLMODEL
         end if
         
!        Describe
!        --------
         call model % Describe
         
      end subroutine InitializeLESModel
!
!/////////////////////////////////////////////////////////////////////////////////////////
!
!           Template procedures
!           -------------------
!
!/////////////////////////////////////////////////////////////////////////////////////////
!
      subroutine LESModel_Initialize(self, controlVariables)
         implicit none
         class(LESModel_t)                     :: self
         class(FTValueDictionary),  intent(in) :: controlVariables

         self % active                = .false.

         !$acc enter data copyin(self)
         !$acc enter data copyin(self % WallModel)
         !$acc enter data copyin(self % LESModel)

      end subroutine LESModel_Initialize

      pure real(kind=RP) function LESModel_ComputeWallEffect (LS,dWall,WallModel)
         !$acc routine seq
         implicit none
         ! class(LESModel_t), intent(in) :: self
         real(kind=RP)    , intent(in) :: LS
         real(kind=RP)    , intent(in) :: dWall
         integer          , intent(in) :: WallModel
         
         select case (WallModel)
            case (LINEAR_WALLMODEL)
               LESModel_ComputeWallEffect = min(LS, dWall * K_VONKARMAN)
            case (NO_WALLMODEL)
               LESModel_ComputeWallEffect = LS
               ! LS is left unmodified if no wall model is selected
         end select
         
      end function LESModel_ComputeWallEffect

      !pure subroutine LESModel_ComputeViscosity (this, delta, dWall, Q, Q_x, Q_y, Q_z, mu)
      !   !$acc routine seq
      !   implicit none
      !   !-arguments---------------------------------------------
      !   class(LESModel_t), intent(in)    :: this
      !   real(kind=RP), intent(in)           :: delta
      !   real(kind=RP), intent(in)           :: dWall
      !   real(kind=RP), intent(in)           :: Q(NCONS)
      !   real(kind=RP), intent(in)           :: Q_x(NGRAD)
      !   real(kind=RP), intent(in)           :: Q_y(NGRAD)
      !   real(kind=RP), intent(in)           :: Q_z(NGRAD)
      !   real(kind=RP), intent(out)          :: mu

      !end subroutine LESModel_ComputeViscosity
      
      subroutine LESModel_Describe(self)
         implicit none
         class(LESModel_t),   intent(in)  :: self


      end subroutine

      subroutine LESModel_Selector(this, delta, dWall, Q, Q_x, Q_y, Q_z, mu)
         !$acc routine seq
         implicit none
         !-arguments---------------------------------------------
         class(LESModel_t), intent(in)       :: this
         real(kind=RP), intent(in)           :: delta
         real(kind=RP), intent(in)           :: dWall
         real(kind=RP), intent(in)           :: Q(NCONS)
         real(kind=RP), intent(in)           :: Q_x(NGRAD)
         real(kind=RP), intent(in)           :: Q_y(NGRAD)
         real(kind=RP), intent(in)           :: Q_z(NGRAD)
         real(kind=RP), intent(out)          :: mu

         select case (this % LESModel)
         case (0) !none
            ! call LESModel_ComputeViscosity(delta, dWall, Q, Q_x, Q_y, Q_z, mu)
         case (1) !Smagorinsky
            call Smagorinsky_ComputeViscosity(delta, dWall, Q, Q_x, Q_y, Q_z, mu, this % C, this % WallModel)
         case (2) !WALE
            call WALE_ComputeViscosity(delta, dWall, Q, Q_x, Q_y, Q_z, mu, this % C)
         case (3) !Vreman
            call Vreman_ComputeViscosity(delta, dWall, Q, Q_x, Q_y, Q_z, mu, this % C)
         end select

      end subroutine LESModel_Selector
!
!//////////////////////////////////////////////////////////////////////////////////////
!
!           Basic Smagorinsky model
!           -----------------------
!
!//////////////////////////////////////////////////////////////////////////////////////
!
      subroutine Smagorinsky_Initialize(self, controlVariables)
         implicit none
         class(Smagorinsky_t)                     :: self
         class(FTValueDictionary),  intent(in) :: controlVariables
!
!        ---------------
!        Local variables
!        ---------------
!
         self % active                = .true.

         if ( controlVariables % containsKey(LESIntensityKey) ) then
            self % C = controlVariables % doublePrecisionValueForKey(LESIntensityKey)

         else
            self % C = 0.2_RP      

         end if

         !$acc enter data copyin(self)
         !$acc enter data copyin(self % C)
         !$acc enter data copyin(self % WallModel)
         !$acc enter data copyin(self % LESModel)

      end subroutine Smagorinsky_Initialize
!
!///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
!
      pure subroutine Smagorinsky_ComputeViscosity (delta, dWall, Q, Q_x, Q_y, Q_z, mu, C, WallModel)
         !$acc routine seq
         implicit none
         !-arguments---------------------------------------------
         ! class(Smagorinsky_t), intent(in)    :: this
         real(kind=RP), intent(in)           :: delta
         real(kind=RP), intent(in)           :: dWall
         real(kind=RP), intent(in)           :: Q(NCONS)
         real(kind=RP), intent(in)           :: Q_x(NGRAD)
         real(kind=RP), intent(in)           :: Q_y(NGRAD)
         real(kind=RP), intent(in)           :: Q_z(NGRAD)
         real(kind=RP), intent(in)           :: C
         integer, intent(in)                 :: WallModel
         real(kind=RP), intent(out)          :: mu
         !-local-variables---------------------------------------
         real(kind=RP)  :: S(NDIM, NDIM)
         real(kind=RP)  :: normS, kappa, LS
         real(kind=RP)  :: U_x(NDIM)
         real(kind=RP)  :: U_y(NDIM)
         real(kind=RP)  :: U_z(NDIM)
         !-------------------------------------------------------
         
         !call getVelocityGradients(Q,Q_x,Q_y,Q_z,U_x,U_y,U_z)
         call getVelocityGradients_State(Q,Q_x,Q_y,Q_z,U_x,U_y,U_z)

!
!        Compute symmetric part of the deformation tensor
!        ------------------------------------------------
         S(:,1) = U_x(1:3)
         S(:,2) = U_y(1:3)
         S(:,3) = U_z(1:3)
         
         S(1,:) = S(1,:) + U_x(1:3)
         S(2,:) = S(2,:) + U_y(1:3)
         S(3,:) = S(3,:) + U_z(1:3)
         
         S = 0.5_RP * S
!
!        Compute the norm of S
!        --------------------- 
         normS = sqrt( 2.0_RP * sum(S*S) )
!
!        Compute viscosity and thermal conductivity
!        ------------------------------------------
         LS = C * delta
         LS = LESModel_ComputeWallEffect(LS,dWall, WallModel)
         mu = Q(IRHO) * POW2(LS) * normS
         
      end subroutine Smagorinsky_ComputeViscosity
!
!///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
!
      subroutine Smagorinsky_Describe(self)
         implicit none
         class(Smagorinsky_t),   intent(in)  :: self

         if ( .not. MPI_Process % isRoot ) return

         write(STD_OUT,*)
         call SubSection_Header("LES Model")
         write(STD_OUT,'(30X,A,A30,A)') "->","LES model: ","Smagorinsky"
         write(STD_OUT,'(30X,A,A30,F10.3)') "->","LES model intensity: ", self % C
         
         select case (self % WallModel)
            case(LINEAR_WALLMODEL)
               write(STD_OUT,'(30X,A,A30,A)') "->","Wall model: ", "linear"
            case(NO_WALLMODEL)
               write(STD_OUT,'(30X,A,A30,A)') "->","Wall model: ", "none"
         end select
         
      end subroutine Smagorinsky_Describe
!
!//////////////////////////////////////////////////////////////////////////////////////
!
!           WALE turbulence model: Wall-Adapting Local Eddy-Viscosity (WALE) Model
!           -----------------------
!            0.55≤Cw≤0.60
!
!//////////////////////////////////////////////////////////////////////////////////////
!
      subroutine WALE_Initialize(self, controlVariables)
         implicit none
         class(WALE_t)                     :: self
         class(FTValueDictionary),  intent(in) :: controlVariables
!
!        ---------------
!        Local variables
!        ---------------
!
         self % active                = .true.

         if ( controlVariables % containsKey(LESIntensityKey) ) then
            self % C = controlVariables % doublePrecisionValueForKey(LESIntensityKey)

         else
            self % C = 0.325_RP      

         end if

         !$acc enter data copyin(self)
         !$acc enter data copyin(self % C)
         !$acc enter data copyin(self % WallModel)
         !$acc enter data copyin(self % LESModel)

      end subroutine WALE_Initialize

      pure subroutine WALE_ComputeViscosity (delta, dWall, Q, Q_x, Q_y, Q_z, mu, C)
         !$acc routine seq
         implicit none
         !-arguments---------------------------------------------
         ! class(WALE_t), intent(in)    :: this
         real(kind=RP), intent(in)           :: delta
         real(kind=RP), intent(in)           :: dWall
         real(kind=RP), intent(in)           :: Q(NCONS)
         real(kind=RP), intent(in)           :: Q_x(NGRAD)
         real(kind=RP), intent(in)           :: Q_y(NGRAD)
         real(kind=RP), intent(in)           :: Q_z(NGRAD)
         real(kind=RP), intent(out)          :: mu
         real(kind=RP), intent(in)           :: C
         !-local-variables---------------------------------------
         real(kind=RP)  :: S(NDIM, NDIM)
         real(kind=RP)  :: gradV2(NDIM, NDIM), gradV(NDIM,NDIM)
         real(kind=RP)  :: Sd(NDIM, NDIM)
         real(kind=RP)  :: normS, normSd, divV2, kappa, LS, rho
         real(kind=RP)  :: U_x(NDIM)
         real(kind=RP)  :: U_y(NDIM)
         real(kind=RP)  :: U_z(NDIM)
         integer        :: i,j
         integer        :: k   ! The third index
         !-------------------------------------------------------
         
         !call getVelocityGradients  (Q,Q_x,Q_y,Q_z,U_x,U_y,U_z)
         call getVelocityGradients_State(Q,Q_x,Q_y,Q_z,U_x,U_y,U_z)

         gradV(1,:) = U_x(1:3)
         gradV(2,:) = U_y(1:3)
         gradV(3,:) = U_z(1:3)
         do i = 1, 3 
            do j = 1, 3 
               S(i,j)      = 0.5_RP*(gradV(i,j)+gradV(j,i))

               gradV2(i,j) = 0
               do k = 1,3
                  gradV2(i,j) = gradV2(i,j) + gradV(i,k)*gradV(k,j) 
               end do

            end do  
         end do  
         divV2 = gradV2(1,1) + gradV2(2,2) + gradV2(3,3)

         normS =  sum(S*S)

         do i = 1, 3 
            do j = 1, 3 
               Sd(i,j) = 0.5_RP*(gradV2(i,j)+gradV2(j,i))
            end do  
         end do  

         Sd(1,1) = Sd(1,1) - 1.0_RP / 3.0_RP * divV2
         Sd(2,2) = Sd(2,2) - 1.0_RP / 3.0_RP * divV2
         Sd(3,3) = Sd(3,3) - 1.0_RP / 3.0_RP * divV2

!//////////////////////////////////////////////////////////////////////////////////////
!        These are the code after expanding the do loop, these codes equal to the above, but hard.
   !   Sd(1,1) = 2.0_RP/3.0*gradV(1,1)**2 + 1.0_RP/3.0*gradV(1,2)*gradV(2,1) + 1.0_RP/3.0*gradV(1,3)*gradV(3,1) - 1.0_RP/3.0*gradV(2,2)**2 - 2.0_RP/3.0*gradV(2,3)*gradV(3,2) - 1.0_RP/3.0*gradV(3,3)**2  
   !   Sd(1,2) = 0.5_RP*gradV(1,1)*gradV(1,2) + 0.5_RP*gradV(1,1)*gradV(2,1) + 0.5_RP*gradV(1,2)*gradV(2,2) + 0.5_RP*gradV(1,3)*gradV(3,2) + 0.5_RP*gradV(2,1)*gradV(2,2) + 0.5_RP*gradV(2,3)*gradV(3,1)  
   !   Sd(1,3) = 0.5_RP*gradV(1,1)*gradV(1,3) + 0.5_RP*gradV(1,1)*gradV(3,1) + 0.5_RP*gradV(1,2)*gradV(2,3) + 0.5_RP*gradV(1,3)*gradV(3,3) + 0.5_RP*gradV(2,1)*gradV(3,2) + 0.5_RP*gradV(3,1)*gradV(3,3)  
   !   Sd(2,1) = 0.5_RP*gradV(1,1)*gradV(1,2) + 0.5_RP*gradV(1,1)*gradV(2,1) + 0.5_RP*gradV(1,2)*gradV(2,2) + 0.5_RP*gradV(1,3)*gradV(3,2) + 0.5_RP*gradV(2,1)*gradV(2,2) + 0.5_RP*gradV(2,3)*gradV(3,1)  
   !   Sd(2,2) = -1.0_RP/3.0*gradV(1,1)**2 + 1.0_RP/3.0*gradV(1,2)*gradV(2,1) - 2.0_RP/3.0*gradV(1,3)*gradV(3,1) + 2.0_RP/3.0*gradV(2,2)**2 + 1.0_RP/3.0*gradV(2,3)*gradV(3,2) - 1.0_RP/3.0*gradV(3,3)**2 
   !   Sd(2,3) = 0.5_RP*gradV(1,2)*gradV(3,1) + 0.5_RP*gradV(1,3)*gradV(2,1) + 0.5_RP*gradV(2,2)*gradV(2,3) + 0.5_RP*gradV(2,2)*gradV(3,2) + 0.5_RP*gradV(2,3)*gradV(3,3) + 0.5_RP*gradV(3,2)*gradV(3,3)  
   !   Sd(3,1) = 0.5_RP*gradV(1,1)*gradV(1,3) + 0.5_RP*gradV(1,1)*gradV(3,1) + 0.5_RP*gradV(1,2)*gradV(2,3) + 0.5_RP*gradV(1,3)*gradV(3,3) + 0.5_RP*gradV(2,1)*gradV(3,2) + 0.5_RP*gradV(3,1)*gradV(3,3)  
   !   Sd(3,2) = 0.5_RP*gradV(1,2)*gradV(3,1) + 0.5_RP*gradV(1,3)*gradV(2,1) + 0.5_RP*gradV(2,2)*gradV(2,3) + 0.5_RP*gradV(2,2)*gradV(3,2) + 0.5_RP*gradV(2,3)*gradV(3,3) + 0.5_RP*gradV(3,2)*gradV(3,3)  
   !   Sd(3,3) = -1.0_RP/3.0*gradV(1,1)**2 - 2.0_RP/3.0*gradV(1,2)*gradV(2,1) + 1.0_RP/3.0*gradV(1,3)*gradV(3,1) - 1.0_RP/3.0*gradV(2,2)**2 + 1.0_RP/3.0*gradV(2,3)*gradV(3,2) + 2.0_RP/3.0*gradV(3,3)**2
!//////////////////////////////////////////////////////////////////////////////////////

         normSd =  sum(Sd*Sd)
         LS = C * delta

         rho = get_rho(Q, dimensionless)

         mu = rho * POW2(LS) * (normSd**(3.0_RP / 2.0_RP) / (normS**(5.0_RP / 2.0_RP)+normSd**(5.0_RP / 4.0_RP)))

         if (normS<1.0e-8_RP .and. normSd<1.0e-8_RP) mu=0.0_RP
         
      end subroutine WALE_ComputeViscosity

      subroutine WALE_Describe(self)
         implicit none
         class(WALE_t),   intent(in)  :: self

         if ( .not. MPI_Process % isRoot ) return

         write(STD_OUT,*)
         call SubSection_Header("LES Model")
         write(STD_OUT,'(30X,A,A30,A)') "->","LES model: ","Wale"
         write(STD_OUT,'(30X,A,A30,F10.3)') "->","LES model intensity: ", self % C
         
         select case (self % WallModel)
            case(NO_WALLMODEL)
               write(STD_OUT,'(30X,A,A30,A)') "->","Wall model: ", "Wale"
            case(LINEAR_WALLMODEL)
               write(STD_OUT,'(30X,A,A30,A)') "->","Wall model: ", "you do not need a linear model with Wale -> deactivate"
         end select
         
      end subroutine WALE_Describe

!
!//////////////////////////////////////////////////////////////////////////////////////
!
!           Vreman
!           -----------------------
!           C = 0.07 (in typical FVM)
!           C = 0.1 in Alya, also recommended by Vreman for highspeed flows 
!//////////////////////////////////////////////////////////////////////////////////////
!
      subroutine Vreman_Initialize(self, controlVariables)
         implicit none
         class(Vreman_t)                     :: self
         class(FTValueDictionary),  intent(in) :: controlVariables
!
!        ---------------
!        Local variables
!        ---------------
!
         self % active                = .true.

         if ( controlVariables % containsKey(LESIntensityKey) ) then
            self % C = controlVariables % doublePrecisionValueForKey(LESIntensityKey)

         else
            self % C = 0.07_RP      

         end if

         !$acc enter data copyin(self)
         !$acc enter data copyin(self % C)
         !$acc enter data copyin(self % WallModel)
         !$acc enter data copyin(self % LESModel)

      end subroutine Vreman_Initialize

      pure subroutine Vreman_ComputeViscosity (delta, dWall, Q, Q_x, Q_y, Q_z, mu, C)
         !$acc routine seq
         implicit none
         !-arguments---------------------------------------------
         ! class(Vreman_t), intent(in)    :: this
         real(kind=RP), intent(in)           :: delta
         real(kind=RP), intent(in)           :: dWall
         real(kind=RP), intent(in)           :: Q(NCONS)
         real(kind=RP), intent(in)           :: Q_x(NGRAD)
         real(kind=RP), intent(in)           :: Q_y(NGRAD)
         real(kind=RP), intent(in)           :: Q_z(NGRAD)
         real(kind=RP), intent(out)          :: mu
         real(kind=RP), intent(in)           :: C
         !-local-variables---------------------------------------
         real(kind=RP)  :: G__ij(NDIM, NDIM)
         real(kind=RP)  :: gradV(NDIM, NDIM)
         real(kind=RP)  :: delta2, alpha, Bbeta, LS
         real(kind=RP)  :: U_x(NDIM)
         real(kind=RP)  :: U_y(NDIM)
         real(kind=RP)  :: U_z(NDIM)
         integer        :: i,j,k
         !-------------------------------------------------------
         
         !call getVelocityGradients  (Q,Q_x,Q_y,Q_z,U_x,U_y,U_z)
         call getVelocityGradients_State(Q,Q_x,Q_y,Q_z,U_x,U_y,U_z)

         delta2 = delta*delta 
         gradV(1,:) = U_x(1:3)
         gradV(2,:) = U_y(1:3)
         gradV(3,:) = U_z(1:3)
         G__ij(:,:) = 0.0_RP

         do i = 1,3
            do j = 1,3
               do k = 1,3
                  G__ij(i,j) = G__ij(i,j) &
                     + (gradV(i,k)*gradV(j,k)*delta2)
               end do
            end do
         end do

         alpha =  sum(gradV*gradV)
         Bbeta = G__ij(1,1) * G__ij(2,2) &
            &  + G__ij(2,2) * G__ij(3,3) &
            &  + G__ij(3,3) * G__ij(1,1) &
            &  - G__ij(1,2) * G__ij(1,2) &
            &  - G__ij(2,3) * G__ij(2,3) &
            &  - G__ij(1,3) * G__ij(1,3)

         if(alpha>1.0e-10_RP) then
            mu = Q(IRHO) * C * sqrt (abs(Bbeta)/alpha)
         else 
            mu = 0.0_RP
         end if
         
      end subroutine Vreman_ComputeViscosity

      subroutine Vreman_Describe(self)
         implicit none
         class(Vreman_t),   intent(in)  :: self

         if ( .not. MPI_Process % isRoot ) return

         write(STD_OUT,*)
         call SubSection_Header("LES Model")
         write(STD_OUT,'(30X,A,A30,A)') "->","LES model: ","Vreman"
         write(STD_OUT,'(30X,A,A30,F10.3)') "->","LES model intensity: ", self % C
         
         select case (self % WallModel)
            case(NO_WALLMODEL)
               write(STD_OUT,'(30X,A,A30,A)') "->","Wall model: ", "Vreman"
            case(LINEAR_WALLMODEL)
               write(STD_OUT,'(30X,A,A30,A)') "->","Wall model: ", "you do not need a linear model with Vreman -> deactivate"
         end select
         
      end subroutine Vreman_Describe

      pure function get_rho(Q, dimensionless_) result(rho)
      !$acc routine seq
         implicit none
         real(kind=RP), intent(in) :: Q(:)              
         type(dimensionless_t),intent(in) :: dimensionless_ 
      
         real(kind=RP) :: rho                          
    
#if defined (NAVIERSTOKES) 
         rho = Q(IRHO)
#elif defined (INCNS)
         rho = Q(INSRHO)
#elif defined (MULTIPHASE)
         rho = dimensionless_%rho(1) * Q(IMC) + dimensionless_%rho(2) * (1.0 - Q(IMC))
! #else
!          print *, "Error: rho computation not valid for physics "
!          stop
#endif
      end function get_rho

end module LESModels
