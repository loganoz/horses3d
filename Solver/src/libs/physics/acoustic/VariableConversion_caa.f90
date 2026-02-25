#include "Includes.h"
module VariableConversion_CAA
   use SMConstants
   use PhysicsStorage_CAA
   use FluidData_CAA
   implicit none

   private
   public   Pressure, PressureDot, PressureBaseFlow_NS, PressureBaseFlow_iNS, PressureBaseFlow_MU
   public   NSGradientVariables_STATE
   ! public   getPrimitiveVariables
   public   getVelocityGradients, getTemperatureGradient, getConservativeGradients
  
   
   interface getConservativeGradients
      module procedure getConservativeGradients_0D, getConservativeGradients_2D, getConservativeGradients_3D
   end interface

   interface getTemperatureGradient
      module procedure getTemperatureGradient_0D, getTemperatureGradient_2D, getTemperatureGradient_3D
   end interface
   
   contains
!
! /////////////////////////////////////////////////////////////////////
!
!@mark -
!---------------------------------------------------------------------
!! Compute the pressure from the state variables
!---------------------------------------------------------------------
!
      PURE function Pressure(Q) RESULT(P)
!
!     ---------
!     Arguments
!     ---------
!
      REAL(KIND=RP), DIMENSION(NCONS), INTENT(IN) :: Q
!
!     ---------------
!     Local Variables
!     ---------------
!
      REAL(KIND=RP) :: P
      
      P = Q(ICAAP)

      end function Pressure
!
! /////////////////////////////////////////////////////////////////////
!---------------------------------------------------------------------
!! Compute the pressure time derivate from the state variables and its time derivatives
!---------------------------------------------------------------------
!
      PURE function PressureDot(Q,QDot) RESULT(PDot)
!
!     ---------
!     Arguments
!     ---------
!
      REAL(KIND=RP), DIMENSION(NCONS), INTENT(IN) :: Q
      REAL(KIND=RP), DIMENSION(NCONS), INTENT(IN) :: QDot
!
!     ---------------
!     Local Variables
!     ---------------
!
      REAL(KIND=RP) :: PDot
      
      PDot = QDot(ICAAP)

      end function PressureDot
!
! /////////////////////////////////////////////////////////////////////
!---------------------------------------------------------------------
!! General procedures for the base flow manipulation
!---------------------------------------------------------------------
! /////////////////////////////////////////////////////////////////////
!
!---------------------------------------------------------------------
!! Compute the pressure from the conservative variables of NS
!---------------------------------------------------------------------
!
      PURE function PressureBaseFlow_NS(Q) RESULT(P)
!
!     ---------
!     Arguments
!     ---------
!
      REAL(KIND=RP), DIMENSION(NCONS), INTENT(IN) :: Q
!
!     ---------------
!     Local Variables
!     ---------------
!
      REAL(KIND=RP) :: P
      
      P = thermodynamics % gammaMinus1*(Q(IRHOE) - 0.5_RP*(Q(IRHOU)**2 + Q(IRHOV)**2 + Q(IRHOW)**2)/Q(IRHO))

      end function PressureBaseFlow_NS
!---------------------------------------------------------------------
!! Compute the pressure from the conservative variables of iNS
!---------------------------------------------------------------------
!
      PURE function PressureBaseFlow_iNS(Q) RESULT(P)
!
!     ---------
!     Arguments
!     ---------
!
      REAL(KIND=RP), DIMENSION(NCONS), INTENT(IN) :: Q
!
!     ---------------
!     Local Variables
!     ---------------
!
      REAL(KIND=RP) :: P
      
      P = Q(INSP)

      end function PressureBaseFlow_iNS
!---------------------------------------------------------------------
!! Compute the pressure from the state variables of MU
!---------------------------------------------------------------------
!
      PURE function PressureBaseFlow_MU(Q) RESULT(P)
!
!     ---------
!     Arguments
!     ---------
!
      REAL(KIND=RP), DIMENSION(NCONS), INTENT(IN) :: Q
!
!     ---------------
!     Local Variables
!     ---------------
!
      REAL(KIND=RP) :: P
      
      P = Q(IMP)

      end function PressureBaseFlow_MU
!
! /////////////////////////////////////////////////////////////////////
!---------------------------------------------------------------------
!! GradientValuesForQ takes the solution (Q) values and returns the
!! quantities of which the gradients of the base flow will be taken.
!---------------------------------------------------------------------
!
      pure subroutine NSGradientVariables_STATE( nEqn, nGrad, Q, U, rho_ )
         implicit none
         integer, intent(in)        :: nEqn, nGrad
         real(kind=RP), intent(in)  :: Q(nEqn)
         real(kind=RP), intent(out) :: U(nGrad)
         real(kind=RP), intent(in), optional :: rho_
!
!        ---------------
!        Local Variables
!        ---------------
!     
         U = Q

      end subroutine NSGradientVariables_STATE
!     
! /////////////////////////////////////////////////////////////////////
!
      pure subroutine getPrimitiveVariablesBaseFlow(U,V)
!
!        **************************************
!           Primitive variables are:
!              V = [rho, u,v,w,p,a2]
!        **************************************
!
         implicit none
         real(kind=RP), intent(in)  :: U(NCONS)
         real(kind=RP), intent(out) :: V(NCONS)
         ! real(kind=RP), intent(out) :: V(NCONS+1)
!
!        ---------------
!        Local variables
!        ---------------
!
         real(kind=RP)  :: invRho

         invRho = 1.0_RP / U(IRHO)

         V(ICAARHO) = U(ICAARHO)
         V(ICAAU) = U(IRHOU) * invRho
         V(ICAAV) = U(IRHOV) * invRho
         V(ICAAW) = U(IRHOW) * invRho
         V(ICAAP) = PressureBaseFlow_NS(U)
         ! V(ICAAA2) = thermodynamics % gamma * V(ICAAP) * invRho

      end subroutine getPrimitiveVariablesBaseFlow
!
!///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
!
!     --------------------------------------
!     Routines to get the velocity gradients
!     --------------------------------------
!
      pure subroutine getVelocityGradients(Q,Q_x,Q_y,Q_z,U_x,U_y,U_z)
         implicit none
         !-arguments---------------------------------------------------
         real(kind=RP), intent(in)  :: Q(NCONS)
         real(kind=RP), intent(in)  :: Q_x(NGRAD), Q_y(NGRAD), Q_z(NGRAD)
         real(kind=RP), intent(out) :: U_x(NDIM), U_y(NDIM), U_z(NDIM)
         !-local-variables---------------------------------------------
         real(kind=RP) :: invRho, invRho2, uDivRho(NDIM)
         !-------------------------------------------------------------

         u_x = Q_x(ICAAU:ICAAW)
         u_y = Q_y(ICAAU:ICAAW)
         u_z = Q_z(ICAAU:ICAAW)

      end subroutine getVelocityGradients
!
!/////////////////////////////////////////////////////////////////////////////
!
!
!///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
!
!     ----------------------------------------------------------
!     Routines to get the temperature gradient
!     -> Currently using the conservative and velocity gradients
!     ----------------------------------------------------------
!
!     ---
!     0D:
!     ---
      pure subroutine getTemperatureGradient_0D(Q,Q_x,Q_y,Q_z,U_x,U_y,U_z,nablaT)
         implicit none
         !-arguments---------------------------------------------------
         real(kind=RP), intent(in)  :: Q(NCONS)
         real(kind=RP), intent(in)  :: Q_x(NGRAD), Q_y(NGRAD), Q_z(NGRAD)
         real(kind=RP), intent(in)  :: U_x(NDIM), U_y(NDIM), U_z(NDIM)
         real(kind=RP), intent(out) :: nablaT(NDIM)
         !-local-variables---------------------------------------------
         real(kind=RP) :: u, v, w, invRho
         !-------------------------------------------------------------
         
         invRho  = 1._RP / Q(IRHO)
         u = Q(IRHOU) / Q(IRHO)
         v = Q(IRHOV) / Q(IRHO)
         w = Q(IRHOW) / Q(IRHO)
         
         nablaT(IX) = thermodynamics % gammaMinus1*dimensionless % gammaM2*(invRho*Q_x(IRHOE) - Q(IRHOE)*invRho*invRho*Q_x(IRHO) - u*u_x(IX)-v*u_x(IY)-w*u_x(IZ))
         nablaT(IY) = thermodynamics % gammaMinus1*dimensionless % gammaM2*(invRho*Q_y(IRHOE) - Q(IRHOE)*invRho*invRho*Q_y(IRHO) - u*u_y(IX)-v*u_y(IY)-w*u_y(IZ))
         nablaT(IZ) = thermodynamics % gammaMinus1*dimensionless % gammaM2*(invRho*Q_z(IRHOE) - Q(IRHOE)*invRho*invRho*Q_z(IRHO) - u*u_z(IX)-v*u_z(IY)-w*u_z(IZ))
      
      end subroutine getTemperatureGradient_0D
!
!/////////////////////////////////////////////////////////////////////////////
!
!     ---
!     2D:
!     ---
      pure subroutine getTemperatureGradient_2D(N,Q,Q_x,Q_y,Q_z,U_x,U_y,U_z,nablaT)
         implicit none
         !-arguments---------------------------------------------------
         integer      , intent(in)  :: N(2)
         real(kind=RP), intent(in)  :: Q  ( NCONS,0:N(1), 0:N(2) )
         real(kind=RP), intent(in)  :: Q_x( NGRAD ,0:N(1), 0:N(2) )
         real(kind=RP), intent(in)  :: Q_y( NGRAD ,0:N(1), 0:N(2) )
         real(kind=RP), intent(in)  :: Q_z( NGRAD ,0:N(1), 0:N(2) )
         real(kind=RP), intent(in)  :: U_x( NDIM ,0:N(1), 0:N(2) )
         real(kind=RP), intent(in)  :: U_y( NDIM ,0:N(1), 0:N(2) )
         real(kind=RP), intent(in)  :: U_z( NDIM ,0:N(1), 0:N(2) )
         real(kind=RP), intent(out) :: nablaT( NDIM ,0:N(1), 0:N(2) )
         !-local-variables---------------------------------------------
         integer :: i,j
         !-------------------------------------------------------------
         
         do j=0, N(2) ; do i=0, N(1)
            call getTemperatureGradient_0D(Q(:,i,j),Q_x(:,i,j),Q_y(:,i,j),Q_z(:,i,j),U_x(:,i,j),U_y(:,i,j),U_z(:,i,j),nablaT(:,i,j))
         end do       ; end do
         
      end subroutine getTemperatureGradient_2D
!
!/////////////////////////////////////////////////////////////////////////////
!
!     ---
!     3D:
!     ---
      pure subroutine getTemperatureGradient_3D(N,Q,Q_x,Q_y,Q_z,U_x,U_y,U_z,nablaT)
         implicit none
         !-arguments---------------------------------------------------
         integer      , intent(in)  :: N(NDIM)
         real(kind=RP), intent(in)  :: Q  ( NCONS ,0:N(1), 0:N(2), 0:N(3) )
         real(kind=RP), intent(in)  :: Q_x( NGRAD ,0:N(1), 0:N(2), 0:N(3) )
         real(kind=RP), intent(in)  :: Q_y( NGRAD ,0:N(1), 0:N(2), 0:N(3) )
         real(kind=RP), intent(in)  :: Q_z( NGRAD ,0:N(1), 0:N(2), 0:N(3) )
         real(kind=RP), intent(out) :: U_x( NDIM  ,0:N(1), 0:N(2), 0:N(3) )
         real(kind=RP), intent(out) :: U_y( NDIM  ,0:N(1), 0:N(2), 0:N(3) )
         real(kind=RP), intent(out) :: U_z( NDIM  ,0:N(1), 0:N(2), 0:N(3) )
         real(kind=RP), intent(out) :: nablaT(NDIM,0:N(1), 0:N(2), 0:N(3) )
         !-local-variables---------------------------------------------
         integer :: i,j,k
         !-------------------------------------------------------------
         
         do k=0, N(3) ; do j=0, N(2) ; do i=0, N(1)
            call getTemperatureGradient_0D(Q(:,i,j,k),Q_x(:,i,j,k),Q_y(:,i,j,k),Q_z(:,i,j,k),U_x(:,i,j,k),U_y(:,i,j,k),U_z(:,i,j,k),nablaT(:,i,j,k))
         end do       ; end do       ; end do
         
      end subroutine getTemperatureGradient_3D
!
!///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
!
!     -----------------------------------------------------------------------
!     Routines to get the conservative gradients from the primitive gradients
!        ( \nabla \rho must already be in Q_x(1), Q_y(1), Q_z(1) )
!     -----------------------------------------------------------------------
!
!     ---
!     0D:
!     ---
      pure subroutine getConservativeGradients_0D(Q,U_x,U_y,U_z,nablaT,Q_x,Q_y,Q_z)
         implicit none
         !-arguments---------------------------------------------------
         real(kind=RP), intent(in)     :: Q(NCONS)
         real(kind=RP), intent(in)     :: U_x(NDIM), U_y(NDIM), U_z(NDIM)
         real(kind=RP), intent(in)     :: nablaT(NDIM)
         real(kind=RP), intent(inout)  :: Q_x(NGRAD), Q_y(NGRAD), Q_z(NGRAD)
         !-local-variables---------------------------------------------
         real(kind=RP) :: u(NDIM), invRho, cons
         !-------------------------------------------------------------
         
         u = Q(IRHOU:IRHOW) / Q(IRHO)
         invRho  = 1._RP / Q(IRHO)
         cons = Q(IRHO) / (thermodynamics % gammaMinus1*dimensionless % gammaM2)
         
         Q_x(IRHOU:IRHOW) = Q(IRHO) * U_x(1:NDIM) + u(1:NDIM) * Q_x(IRHO)
         Q_y(IRHOU:IRHOW) = Q(IRHO) * U_y(1:NDIM) + u(1:NDIM) * Q_y(IRHO)
         Q_z(IRHOU:IRHOW) = Q(IRHO) * U_z(1:NDIM) + u(1:NDIM) * Q_z(IRHO)
         
         Q_x(IRHOE) = cons * nablaT(IX) + Q(IRHOE) * invRho * Q_x(IRHO) + Q(IRHOU) * U_x(IX) + Q(IRHOV) * U_x(IY) + Q(IRHOV) * U_x(IZ)
         Q_y(IRHOE) = cons * nablaT(IY) + Q(IRHOE) * invRho * Q_y(IRHO) + Q(IRHOU) * U_y(IX) + Q(IRHOV) * U_y(IY) + Q(IRHOV) * U_y(IZ)
         Q_z(IRHOE) = cons * nablaT(IZ) + Q(IRHOE) * invRho * Q_z(IRHO) + Q(IRHOU) * U_z(IX) + Q(IRHOV) * U_z(IY) + Q(IRHOV) * U_z(IZ)
         
      end subroutine getConservativeGradients_0D
!     ---
!     2D:
!     ---
      pure subroutine getConservativeGradients_2D(N,Q,U_x,U_y,U_z,nablaT,Q_x,Q_y,Q_z)
         implicit none
         !-arguments---------------------------------------------------
         integer      , intent(in)     :: N     (2)
         real(kind=RP), intent(in)     :: Q     (NCONS, 0:N(1), 0:N(2))
         real(kind=RP), intent(in)     :: U_x   (NDIM , 0:N(1), 0:N(2))
         real(kind=RP), intent(in)     :: U_y   (NDIM , 0:N(1), 0:N(2))
         real(kind=RP), intent(in)     :: U_z   (NDIM , 0:N(1), 0:N(2))
         real(kind=RP), intent(in)     :: nablaT(NDIM , 0:N(1), 0:N(2))
         real(kind=RP), intent(inout)  :: Q_x   (NGRAD, 0:N(1), 0:N(2))
         real(kind=RP), intent(inout)  :: Q_y   (NGRAD, 0:N(1), 0:N(2))
         real(kind=RP), intent(inout)  :: Q_z   (NGRAD, 0:N(1), 0:N(2))
         !-local-variables---------------------------------------------
         integer :: i,j
         !-------------------------------------------------------------
         
         do j=0, N(2) ; do i=0, N(1)
            call getConservativeGradients_0D(Q(:,i,j),U_x(:,i,j),U_y(:,i,j),U_z(:,i,j),nablaT(:,i,j),Q_x(:,i,j),Q_y(:,i,j),Q_z(:,i,j))
         end do       ; end do
         
      end subroutine getConservativeGradients_2D
!     ---
!     3D:
!     ---
      pure subroutine getConservativeGradients_3D(N,Q,U_x,U_y,U_z,nablaT,Q_x,Q_y,Q_z)
         implicit none
         !-arguments---------------------------------------------------
         integer      , intent(in)     :: N     (3)
         real(kind=RP), intent(in)     :: Q     (NCONS, 0:N(1), 0:N(2), 0:N(3))
         real(kind=RP), intent(in)     :: U_x   (NDIM , 0:N(1), 0:N(2), 0:N(3))
         real(kind=RP), intent(in)     :: U_y   (NDIM , 0:N(1), 0:N(2), 0:N(3))
         real(kind=RP), intent(in)     :: U_z   (NDIM , 0:N(1), 0:N(2), 0:N(3))
         real(kind=RP), intent(in)     :: nablaT(NDIM , 0:N(1), 0:N(2), 0:N(3))
         real(kind=RP), intent(inout)  :: Q_x   (NGRAD, 0:N(1), 0:N(2), 0:N(3))
         real(kind=RP), intent(inout)  :: Q_y   (NGRAD, 0:N(1), 0:N(2), 0:N(3))
         real(kind=RP), intent(inout)  :: Q_z   (NGRAD, 0:N(1), 0:N(2), 0:N(3))
         !-local-variables---------------------------------------------
         integer :: i,j, k
         !-------------------------------------------------------------
         
         do k=0, N(3) ; do j=0, N(2) ; do i=0, N(1)
            call getConservativeGradients_0D(Q(:,i,j,k),U_x(:,i,j,k),U_y(:,i,j,k),U_z(:,i,j,k),nablaT(:,i,j,k),Q_x(:,i,j,k),Q_y(:,i,j,k),Q_z(:,i,j,k))
         end do       ; end do       ; end do
         
      end subroutine getConservativeGradients_3D



      ! pure subroutine computeSoundVelocity_NS()
      !       implicit none


      ! end subroutine computeSoundVelocity_NS

end module VariableConversion_CAA
