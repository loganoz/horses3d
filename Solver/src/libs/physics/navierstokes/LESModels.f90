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
   use PhysicsStorage_NS
   use FTValueDictionaryClass
   use Physics_NSKeywordsModule
   use MPI_Process_Info
   use Headers
   use Utilities                 , only: toLower
   use FluidData_NS
   use VariableConversion_NS     , only: getVelocityGradients, getTemperatureGradient
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

            case default
               write(STD_OUT,'(A,A,A)') "LES Model ",trim(modelName), " is not implemented."
               print*, "Available options are:"
               print*, "   * None (default)"
               print*, "   * Smagorinsky"
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

            case default
               write(STD_OUT,'(A,A,A)') "Wall model ",trim(modelName), " is not implemented."
               print*, "Available options are:"
               print*, "   * Linear (default)"
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

      pure subroutine LESModel_ComputeViscosity (this, delta, dWall, Q, Q_x, Q_y, Q_z, muSmag)
         implicit none
         !-arguments---------------------------------------------
         class(LESModel_t), intent(in)    :: this
         real(kind=RP), intent(in)           :: delta
         real(kind=RP), intent(in)           :: dWall
         real(kind=RP), intent(in)           :: Q(NCONS)
         real(kind=RP), intent(in)           :: Q_x(NGRAD)
         real(kind=RP), intent(in)           :: Q_y(NGRAD)
         real(kind=RP), intent(in)           :: Q_z(NGRAD)
         real(kind=RP), intent(out)          :: muSmag

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
      pure subroutine Smagorinsky_ComputeViscosity (this, delta, dWall, Q, Q_x, Q_y, Q_z, muSmag)
         implicit none
         !-arguments---------------------------------------------
         class(Smagorinsky_t), intent(in)    :: this
         real(kind=RP), intent(in)           :: delta
         real(kind=RP), intent(in)           :: dWall
         real(kind=RP), intent(in)           :: Q(NCONS)
         real(kind=RP), intent(in)           :: Q_x(NGRAD)
         real(kind=RP), intent(in)           :: Q_y(NGRAD)
         real(kind=RP), intent(in)           :: Q_z(NGRAD)
         real(kind=RP), intent(out)          :: muSmag
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
         muSmag = Q(IRHO) * POW2(LS) * normS
         
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

end module LESModels
