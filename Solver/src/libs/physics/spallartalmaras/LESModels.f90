!
!//////////////////////////////////////////////////////
!
!   @File:    LESModels.f90
!   @Author:  Juan Manzanero (juan.manzanero@upm.es)
!   @Created: Sun Jan 14 13:23:10 2018
!   @Last revision date: Thu May  2 09:41:43 2019
!   @Last revision author: Andrés Rueda (am.rueda@upm.es)
!   @Last revision commit: 67c9993eab2425db318bd6a45ef48d4abba673b7
!
!//////////////////////////////////////////////////////
!
#include "Includes.h"
module LESModels
   use SMConstants
   use PhysicsStorage_NSSA
   use FTValueDictionaryClass
   use Physics_NSSAKeywordsModule
   use MPI_Process_Info
   use Headers
   use Utilities                 , only: toLower
   use FluidData_NSSA
   use VariableConversion_NSSA     , only: getVelocityGradients, getTemperatureGradient
   implicit none

   private
   public LESModel, InitializeLESModel, Smagorinsky_t
   
!  Keywords
!  --------
   character(len=*), parameter   :: LESIntensityKey = "les model intensity"
   character(len=*), parameter   :: WallModelKey    = "wall model"
   
!  Model parameters
!  ----------------
   real(kind=RP)   , parameter   :: K_VONKARMAN     = 0.4_RP
   
!  Wall models
!  -----------
   integer         , parameter   :: NO_WALLMODEL       = 0
   integer         , parameter   :: LINEAR_WALLMODEL   = 1
   
   type LESModel_t
      logical  :: active
      logical  :: requiresWallDistances = .FALSE.
      integer  :: WallModel
      contains
         procedure            :: Initialize         => LESModel_Initialize
         procedure, private   :: ComputeWallEffect  => LESModel_ComputeWallEffect
         procedure            :: ComputeViscosity   => LESModel_ComputeViscosity
         procedure            :: Describe           => LESModel_Describe
   end type LESModel_t

   type, extends(LESModel_t)  :: Smagorinsky_t
      real(kind=RP)  :: CS
      contains
         procedure          :: Initialize         => Smagorinsky_Initialize
         procedure          :: Describe           => Smagorinsky_Describe
         procedure          :: ComputeViscosity   => Smagorinsky_ComputeViscosity
   end type Smagorinsky_t

   type, extends(LESModel_t)  :: WALE_t
      real(kind=RP)  :: Cw
      contains
         procedure          :: Initialize         => WALE_Initialize
         procedure          :: Describe           => WALE_Describe
         procedure          :: ComputeViscosity   => WALE_ComputeViscosity
   end type WALE_t


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

            case ("smagorinsky")
               if (.not. allocated(model)) allocate(Smagorinsky_t  :: model)

            case ("wale")
               if (.not. allocated(model)) allocate(WALE_t  :: model)

            case default
               write(STD_OUT,'(A,A,A)') "LES Model ",trim(modelName), " is not implemented."
               print*, "Available options are:"
               print*, "   * None (default)"
               print*, "   * Smagorinsky"
               print*, "   * Wale"
               errorMessage(STD_OUT)
               stop

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

            case ("linear")
               model % WallModel             = LINEAR_WALLMODEL
               model % requiresWallDistances = .true.

            case ("Wale")
               model % WallModel = NO_WALLMODEL
               model % requiresWallDistances = .true.

            case default
               write(STD_OUT,'(A,A,A)') "Wall model ",trim(modelName), " is not implemented."
               print*, "Available options are:"
               print*, "   * Linear (default)"
               print*, "   * Wale"
               print*, "   * None"
               errorMessage(STD_OUT)
               stop

            end select
            
         else
            model % WallModel = LINEAR_WALLMODEL
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

      end subroutine LESModel_Initialize

      pure real(kind=RP) function LESModel_ComputeWallEffect (self,LS,dWall)
         implicit none
         class(LESModel_t), intent(in) :: self
         real(kind=RP)    , intent(in) :: LS
         real(kind=RP)    , intent(in) :: dWall
         
         select case (self % WallModel)
            case (LINEAR_WALLMODEL)
               LESModel_ComputeWallEffect = min(LS, dWall * K_VONKARMAN)
         end select
         
      end function LESModel_ComputeWallEffect

      pure subroutine LESModel_ComputeViscosity (this, delta, dWall, Q, Q_x, Q_y, Q_z, mu)
         implicit none
         !-arguments---------------------------------------------
         class(LESModel_t), intent(in)    :: this
         real(kind=RP), intent(in)           :: delta
         real(kind=RP), intent(in)           :: dWall
         real(kind=RP), intent(in)           :: Q(NCONS)
         real(kind=RP), intent(in)           :: Q_x(NGRAD)
         real(kind=RP), intent(in)           :: Q_y(NGRAD)
         real(kind=RP), intent(in)           :: Q_z(NGRAD)
         real(kind=RP), intent(out)          :: mu

      end subroutine LESModel_ComputeViscosity
      
      subroutine LESModel_Describe(self)
         implicit none
         class(LESModel_t),   intent(in)  :: self


      end subroutine
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
            self % CS = controlVariables % doublePrecisionValueForKey(LESIntensityKey)

         else
            self % CS = 0.2_RP      

         end if

      end subroutine Smagorinsky_Initialize
!
!///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
!
      pure subroutine Smagorinsky_ComputeViscosity (this, delta, dWall, Q, Q_x, Q_y, Q_z, mu)
         implicit none
         !-arguments---------------------------------------------
         class(Smagorinsky_t), intent(in)    :: this
         real(kind=RP), intent(in)           :: delta
         real(kind=RP), intent(in)           :: dWall
         real(kind=RP), intent(in)           :: Q(NCONS)
         real(kind=RP), intent(in)           :: Q_x(NGRAD)
         real(kind=RP), intent(in)           :: Q_y(NGRAD)
         real(kind=RP), intent(in)           :: Q_z(NGRAD)
         real(kind=RP), intent(out)          :: mu
         !-local-variables---------------------------------------
         real(kind=RP)  :: S(NDIM, NDIM)
         real(kind=RP)  :: normS, kappa, LS
         real(kind=RP)  :: U_x(NDIM)
         real(kind=RP)  :: U_y(NDIM)
         real(kind=RP)  :: U_z(NDIM)
         !-------------------------------------------------------
         
         call getVelocityGradients(Q,Q_x,Q_y,Q_z,U_x,U_y,U_z)
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
         LS = this % CS * delta
         LS = this % ComputeWallEffect(LS,dWall)
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
         write(STD_OUT,'(30X,A,A30,F10.3)') "->","LES model intensity: ", self % CS
         
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
            self % Cw = controlVariables % doublePrecisionValueForKey(LESIntensityKey)

         else
            self % Cw = 0.325_RP      

         end if

      end subroutine WALE_Initialize

      pure subroutine WALE_ComputeViscosity (this, delta, dWall, Q, Q_x, Q_y, Q_z, mu)
         implicit none
         !-arguments---------------------------------------------
         class(WALE_t), intent(in)    :: this
         real(kind=RP), intent(in)           :: delta
         real(kind=RP), intent(in)           :: dWall
         real(kind=RP), intent(in)           :: Q(NCONS)
         real(kind=RP), intent(in)           :: Q_x(NGRAD)
         real(kind=RP), intent(in)           :: Q_y(NGRAD)
         real(kind=RP), intent(in)           :: Q_z(NGRAD)
         real(kind=RP), intent(out)          :: mu
         !-local-variables---------------------------------------
         real(kind=RP)  :: S(NDIM, NDIM)
         real(kind=RP)  :: Sd(NDIM, NDIM)
         real(kind=RP)  :: normS, normSd, divV, kappa, LS
         real(kind=RP)  :: U_x(NDIM)
         real(kind=RP)  :: U_y(NDIM)
         real(kind=RP)  :: U_z(NDIM)
         integer        :: m
         !-------------------------------------------------------
         
         call getVelocityGradients  (Q,Q_x,Q_y,Q_z,U_x,U_y,U_z)
         
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

         divV = S(1,1) + S(2,2) + S(3,3)

!        Compute the norm of S
!        --------------------- 
         normS =  sum(S*S)

!        Remove the volumetric deformation tensor with squared gradients
!        ----------------------------------------
	 do m = 1, 3 
         Sd(m,1) = POW2(U_x(m))
         Sd(m,2) = POW2(U_y(m))
         Sd(m,3) = POW2(U_z(m))
	 end do  

	 do m = 1, 3
         Sd(1,m) = Sd(1,m) + POW2(U_x(m))
         Sd(2,m) = Sd(2,m) + POW2(U_y(m))
         Sd(3,m) = Sd(3,m) + POW2(U_z(m))
	 end do 

         Sd = 0.5_RP * Sd

	      Sd(1,1) = Sd(1,1) - 1.0_RP / 3.0_RP * POW2(divV)
         Sd(2,2) = Sd(2,2) - 1.0_RP / 3.0_RP * POW2(divV)
         Sd(3,3) = Sd(3,3) - 1.0_RP / 3.0_RP * POW2(divV)

!        Compute the norm of Sd
!        --------------------- 
         normSd =  sum(Sd*Sd)
!
!        Compute viscosity and thermal conductivity
!        ------------------------------------------
         LS = min(this % Cw * delta, dWall * K_VONKARMAN)
         LS = this % ComputeWallEffect(LS,dWall)
         
         mu = Q(IRHO) * POW2(LS) * (normSd**(3.0_RP / 2.0_RP) / (normS**(5.0_RP / 2.0_RP)+normSd**(5.0_RP / 4.0_RP)))
         
      end subroutine WALE_ComputeViscosity

      subroutine WALE_Describe(self)
         implicit none
         class(WALE_t),   intent(in)  :: self

         if ( .not. MPI_Process % isRoot ) return

         write(STD_OUT,*)
         call SubSection_Header("LES Model")
         write(STD_OUT,'(30X,A,A30,A)') "->","LES model: ","Wale"
         write(STD_OUT,'(30X,A,A30,F10.3)') "->","LES model intensity: ", self % Cw
         
         select case (self % WallModel)
            case(NO_WALLMODEL)
               write(STD_OUT,'(30X,A,A30,A)') "->","Wall model: ", "Wale"
            case(LINEAR_WALLMODEL)
               write(STD_OUT,'(30X,A,A30,A)') "->","Wall model: ", "you do not need a linear model wiht Wale -> deactivate"
         end select
         
      end subroutine WALE_Describe

end module LESModels
