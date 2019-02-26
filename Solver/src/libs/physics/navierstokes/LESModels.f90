!
!//////////////////////////////////////////////////////
!
!   @File:    LESModels.f90
!   @Author:  Juan Manzanero (juan.manzanero@upm.es)
!   @Created: Sun Jan 14 13:23:10 2018
!   @Last revision date: Tue Feb 26 18:17:38 2019
!   @Last revision author: Andrés Rueda (am.rueda@upm.es)
!   @Last revision commit: 2bf19fe8b635bd13e65827bd1023dd656a765b42
!
!//////////////////////////////////////////////////////
!
!
!//////////////////////////////////////////////////////
!
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

   real(kind=RP), parameter      :: K_VONKARMAN = 0.4_RP
   character(len=*), parameter   :: LESIntensityKey = "les model intensity"

   type LESModel_t
      logical  :: active
      logical  :: requiresWallDistances
      contains
         procedure            :: Initialize         => LESModel_Initialize
         generic              :: ComputeSGSTensor   => ComputeSGSTensor0D, ComputeSGSTensor2D, ComputeSGSTensor3D
         procedure, private   :: ComputeSGSTensor3D => LESModel_ComputeSGSTensor3D
         procedure, private   :: ComputeSGSTensor2D => LESModel_ComputeSGSTensor2D
         procedure, private   :: ComputeSGSTensor0D => LESModel_ComputeSGSTensor0D
         procedure            :: Describe           => LESModel_Describe
   end type LESModel_t

   type, extends(LESModel_t)  :: Smagorinsky_t
      real(kind=RP)  :: CS
      contains
         procedure          :: Initialize         => Smagorinsky_Initialize
         procedure, private :: ComputeSGSTensor3D => Smagorinsky_ComputeSGSTensor3D
         procedure, private :: ComputeSGSTensor2D => Smagorinsky_ComputeSGSTensor2D
         procedure, private :: ComputeSGSTensor0D => Smagorinsky_ComputeSGSTensor0D
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
         self % requiresWallDistances = .false.

      end subroutine LESModel_Initialize

      subroutine LESModel_ComputeSGSTensor3D(self, delta, N, dWall, Q, Q_x, Q_y, Q_z, tau, qSGS)
         implicit none
         class(LESModel_t), intent(in) :: self
         real(kind=RP), intent(in)     :: delta
         integer,       intent(in)     :: N(3)
         real(kind=RP), intent(in)     :: dWall(0:N(1), 0:N(2), 0:N(3))
         real(kind=RP), intent(in)     :: Q  (NCONS, 0:N(1), 0:N(2), 0:N(3))
         real(kind=RP), intent(in)     :: Q_x(NGRAD, 0:N(1), 0:N(2), 0:N(3))
         real(kind=RP), intent(in)     :: Q_y(NGRAD, 0:N(1), 0:N(2), 0:N(3))
         real(kind=RP), intent(in)     :: Q_z(NGRAD, 0:N(1), 0:N(2), 0:N(3))
         real(kind=RP), intent(out)    :: tau(NDIM, NDIM, 0:N(1), 0:N(2), 0:N(3))
         real(kind=RP), intent(out)    :: qSGS(NDIM, 0:N(1), 0:N(2), 0:N(3))
           
         tau = 0.0_RP
         qSGS= 0.0_RP

      end subroutine LESModel_ComputeSGSTensor3D

      subroutine LESModel_ComputeSGSTensor2D(self, delta, N, dWall, Q, Q_x, Q_y, Q_z, tau, qSGS)
         implicit none
         class(LESModel_t), intent(in) :: self
         real(kind=RP), intent(in)     :: delta
         integer,       intent(in)     :: N(2)
         real(kind=RP), intent(in)     :: dWall(0:N(1), 0:N(2))
         real(kind=RP), intent(in)     :: Q(NCONS, 0:N(1), 0:N(2))
         real(kind=RP), intent(in)     :: Q_x(NGRAD, 0:N(1), 0:N(2))
         real(kind=RP), intent(in)     :: Q_y(NGRAD, 0:N(1), 0:N(2))
         real(kind=RP), intent(in)     :: Q_z(NGRAD, 0:N(1), 0:N(2))
         real(kind=RP), intent(out)    :: tau(NDIM, NDIM, 0:N(1), 0:N(2))
         real(kind=RP), intent(out)    :: qSGS(NDIM, 0:N(1), 0:N(2))
           
         tau = 0.0_RP
         qSGS= 0.0_RP

      end subroutine LESModel_ComputeSGSTensor2D

      subroutine LESModel_ComputeSGSTensor0D(self, delta, dWall, Q, Q_x, Q_y, Q_z, tau, qSGS)
         implicit none
         class(LESModel_t), intent(in) :: self
         real(kind=RP), intent(in)     :: delta
         real(kind=RP), intent(in)     :: dWall
         real(kind=RP), intent(in)     :: Q(NCONS)
         real(kind=RP), intent(in)     :: Q_x(NGRAD)
         real(kind=RP), intent(in)     :: Q_y(NGRAD)
         real(kind=RP), intent(in)     :: Q_z(NGRAD)
         real(kind=RP), intent(out)    :: tau(NDIM, NDIM)
         real(kind=RP), intent(out)    :: qSGS(NDIM)

         tau = 0.0_RP
         qSGS= 0.0_RP

      end subroutine LESModel_ComputeSGSTensor0D

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
         self % requiresWallDistances = .true.

         if ( controlVariables % containsKey(LESIntensityKey) ) then
            self % CS = controlVariables % doublePrecisionValueForKey(LESIntensityKey)

         else
            self % CS = 0.2_RP      

         end if

      end subroutine Smagorinsky_Initialize

      subroutine Smagorinsky_ComputeSGSTensor3D(self, delta, N, dWall, Q, Q_x, Q_y, Q_z, tau, qSGS)
         implicit none
         class(Smagorinsky_t), intent(in) :: self
         real(kind=RP), intent(in)     :: delta
         integer,       intent(in)     :: N(3)
         real(kind=RP), intent(in)     :: dWall(0:N(1), 0:N(2), 0:N(3))
         real(kind=RP), intent(in)     :: Q  (NCONS, 0:N(1), 0:N(2), 0:N(3))
         real(kind=RP), intent(in)     :: Q_x(NGRAD, 0:N(1), 0:N(2), 0:N(3))
         real(kind=RP), intent(in)     :: Q_y(NGRAD, 0:N(1), 0:N(2), 0:N(3))
         real(kind=RP), intent(in)     :: Q_z(NGRAD, 0:N(1), 0:N(2), 0:N(3))
         real(kind=RP), intent(out)    :: tau(NDIM, NDIM, 0:N(1), 0:N(2), 0:N(3))
         real(kind=RP), intent(out)    :: qSGS(NDIM, 0:N(1), 0:N(2), 0:N(3))
!
!        ---------------
!        Local variables
!        ---------------
!
         integer     :: i, j, k
         real(kind=RP)  :: S(NDIM, NDIM)
         real(kind=RP)  :: normS, divV, mu, kappa, LS
         real(kind=RP) :: U_x(NDIM, 0:N(1), 0:N(2), 0:N(3))
         real(kind=RP) :: U_y(NDIM, 0:N(1), 0:N(2), 0:N(3))
         real(kind=RP) :: U_z(NDIM, 0:N(1), 0:N(2), 0:N(3))
         real(kind=RP) :: nablaT(NDIM, 0:N(1), 0:N(2), 0:N(3))
         
         call getVelocityGradients  (N,Q,Q_x,Q_y,Q_z,U_x,U_y,U_z)
         call getTemperatureGradient(N,Q,Q_x,Q_y,Q_z,U_x,U_y,U_z,nablaT)

         do k = 0, N(3) ; do j = 0, N(2)  ; do i = 0, N(1)
!
!           Compute symmetric part of the deformation tensor
!           ------------------------------------------------
            S(:,1) = U_x(1:3, i, j, k)
            S(:,2) = U_y(1:3, i, j, k)
            S(:,3) = U_z(1:3, i, j, k)

            S(1,:) = S(1,:) + U_x(1:3,i,j,k)
            S(2,:) = S(2,:) + U_y(1:3,i,j,k)
            S(3,:) = S(3,:) + U_z(1:3,i,j,k)

            S = 0.5_RP * S

            divV = S(1,1) + S(2,2) + S(3,3)
!
!           Compute the norm of S
!           --------------------- 
            normS = sqrt( 2.0_RP * sum(S*S) )
!
!           Compute viscosity and thermal conductivity
!           ------------------------------------------
            LS = min(self % CS * delta, dWall(i,j,k) * K_VONKARMAN)
            mu = POW2(LS) * normS
            kappa = mu / (thermodynamics % gammaMinus1 * POW2(dimensionless % Mach) * dimensionless % Pr)
!
!           Remove the volumetric deformation tensor
!           ----------------------------------------
            S(1,1) = S(1,1) - 1.0_RP / 3.0_RP * divV
            S(2,2) = S(2,2) - 1.0_RP / 3.0_RP * divV
            S(3,3) = S(3,3) - 1.0_RP / 3.0_RP * divV
!
!           Compute the SGS tensor and heat flux
!           ------------------------------------
            tau(:,:,i,j,k) = -2.0_RP * mu * S

            qSGS(1,i,j,k) = -kappa * nablaT(IX,i,j,k)
            qSGS(2,i,j,k) = -kappa * nablaT(IY,i,j,k)
            qSGS(3,i,j,k) = -kappa * nablaT(IZ,i,j,k)

         end do         ; end do          ; end do
         
      end subroutine Smagorinsky_ComputeSGSTensor3D

      subroutine Smagorinsky_ComputeSGSTensor2D(self, delta, N, dWall, Q, Q_x, Q_y, Q_z, tau, qSGS)
         implicit none
         class(Smagorinsky_t), intent(in) :: self
         real(kind=RP), intent(in)     :: delta
         integer,       intent(in)     :: N(2)
         real(kind=RP), intent(in)     :: dWall(0:N(1), 0:N(2))
         real(kind=RP), intent(in)     :: Q(NCONS, 0:N(1), 0:N(2))
         real(kind=RP), intent(in)     :: Q_x(NGRAD, 0:N(1), 0:N(2))
         real(kind=RP), intent(in)     :: Q_y(NGRAD, 0:N(1), 0:N(2))
         real(kind=RP), intent(in)     :: Q_z(NGRAD, 0:N(1), 0:N(2))
         real(kind=RP), intent(out)    :: tau(NDIM, NDIM, 0:N(1), 0:N(2))
         real(kind=RP), intent(out)    :: qSGS(NDIM, 0:N(1), 0:N(2))
!
!        ---------------
!        Local variables
!        ---------------
!
         integer       :: i, j
         real(kind=RP) :: S(NDIM, NDIM)
         real(kind=RP) :: normS, divV, mu, kappa, LS
         real(kind=RP) :: U_x(NDIM, 0:N(1), 0:N(2))
         real(kind=RP) :: U_y(NDIM, 0:N(1), 0:N(2))
         real(kind=RP) :: U_z(NDIM, 0:N(1), 0:N(2))
         real(kind=RP) :: nablaT(NDIM, 0:N(1), 0:N(2))
         
         call getVelocityGradients  (N,Q,Q_x,Q_y,Q_z,U_x,U_y,U_z)
         call getTemperatureGradient(N,Q,Q_x,Q_y,Q_z,U_x,U_y,U_z,nablaT)
         
         do j = 0, N(2)  ; do i = 0, N(1)
!
!           Compute symmetric part of the deformation tensor
!           ------------------------------------------------
            S(:,1) = U_x(1:3, i, j)
            S(:,2) = U_y(1:3, i, j)
            S(:,3) = U_z(1:3, i, j)

            S(1,:) = S(1,:) + U_x(1:3,i,j)
            S(2,:) = S(2,:) + U_y(1:3,i,j)
            S(3,:) = S(3,:) + U_z(1:3,i,j)

            S = 0.5_RP * S

            divV = S(1,1) + S(2,2) + S(3,3)
!
!           Compute the norm of S
!           --------------------- 
            normS = sqrt( 2.0_RP * sum(S*S) )
!
!           Compute viscosity and thermal conductivity
!           ------------------------------------------
            LS = min(self % CS * delta, dWall(i,j) * K_VONKARMAN)
            mu = POW2(LS) * normS
            kappa = mu / (thermodynamics % gammaMinus1 * POW2(dimensionless % Mach) * dimensionless % Pr)
!
!           Remove the volumetric deformation tensor
!           ----------------------------------------
            S(1,1) = S(1,1) - 1.0_RP / 3.0_RP * divV
            S(2,2) = S(2,2) - 1.0_RP / 3.0_RP * divV
            S(3,3) = S(3,3) - 1.0_RP / 3.0_RP * divV
!
!           Compute the SGS tensor and heat flux
!           ------------------------------------
            tau(:,:,i,j) = -2.0_RP * mu * S

            qSGS(1,i,j) = -kappa * nablaT(IX,i,j)
            qSGS(2,i,j) = -kappa * nablaT(IY,i,j)
            qSGS(3,i,j) = -kappa * nablaT(IZ,i,j)

         end do          ; end do
         
      end subroutine Smagorinsky_ComputeSGSTensor2D

      subroutine Smagorinsky_ComputeSGSTensor0D(self, delta, dWall, Q, Q_x, Q_y, Q_z, tau, qSGS)
         implicit none
         class(Smagorinsky_t), intent(in) :: self
         real(kind=RP), intent(in)     :: delta
         real(kind=RP), intent(in)     :: dWall
         real(kind=RP), intent(in)     :: Q(NCONS)
         real(kind=RP), intent(in)     :: Q_x(NGRAD)
         real(kind=RP), intent(in)     :: Q_y(NGRAD)
         real(kind=RP), intent(in)     :: Q_z(NGRAD)
         real(kind=RP), intent(out)    :: tau(NDIM, NDIM)
         real(kind=RP), intent(out)    :: qSGS(NDIM)
!
!        ---------------
!        Local variables
!        ---------------
!
         real(kind=RP)  :: S(NDIM, NDIM)
         real(kind=RP)  :: normS, divV, mu, kappa, LS
         real(kind=RP) :: U_x(NDIM)
         real(kind=RP) :: U_y(NDIM)
         real(kind=RP) :: U_z(NDIM)
         real(kind=RP) :: nablaT(NDIM)
         
         call getVelocityGradients  (Q,Q_x,Q_y,Q_z,U_x,U_y,U_z)
         call getTemperatureGradient(Q,Q_x,Q_y,Q_z,U_x,U_y,U_z,nablaT)
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
!
!        Compute the norm of S
!        --------------------- 
         normS = sqrt( 2.0_RP * sum(S*S) )
!
!        Compute viscosity and thermal conductivity
!        ------------------------------------------
         LS = min(self % CS * delta, dWall * K_VONKARMAN)
         mu = POW2(LS) * normS
         kappa = mu / (thermodynamics % gammaMinus1 * POW2(dimensionless % Mach) * dimensionless % Pr)
!
!        Remove the volumetric deformation tensor
!        ----------------------------------------
         S(1,1) = S(1,1) - 1.0_RP / 3.0_RP * divV
         S(2,2) = S(2,2) - 1.0_RP / 3.0_RP * divV
         S(3,3) = S(3,3) - 1.0_RP / 3.0_RP * divV
!
!        Compute the SGS tensor
!        ----------------------
         tau = -2.0_RP * mu * S

         qSGS(1) = -kappa * nablaT(IX)
         qSGS(2) = -kappa * nablaT(IY)
         qSGS(3) = -kappa * nablaT(IZ)

      end subroutine Smagorinsky_ComputeSGSTensor0D
!
!///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
!
      subroutine Smagorinsky_ComputeViscosity (this, delta, dWall, Q, Q_x, Q_y, Q_z, muSmag)
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
!
!        Compute the norm of S
!        --------------------- 
         normS = sqrt( 2.0_RP * sum(S*S) )
!
!        Compute viscosity and thermal conductivity
!        ------------------------------------------
         LS = min(this % CS * delta, dWall * K_VONKARMAN)
         muSmag = POW2(LS) * normS
         
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

      end subroutine Smagorinsky_Describe

end module LESModels
