!
!//////////////////////////////////////////////////////
!
!   @File:    Physics_iNS.f90
!   @Author:  Juan Manzanero (juan.manzanero@upm.es)
!   @Created: Tue Jun 19 17:39:26 2018
!   @Last revision date: Wed Aug  1 15:48:18 2018
!   @Last revision author: Juan Manzanero (juan.manzanero@upm.es)
!   @Last revision commit: f358d5850cf9ae49fb85272ef0ea077425d7ed8b
!
!//////////////////////////////////////////////////////
!
#include "Includes.h"
!  **************
   module Physics_iNS
!  **************
!
      use SMConstants
      use PhysicsStorage_iNS
      use VariableConversion_iNS
      use FluidData_iNS
      implicit none

      private
      public  iEulerFlux, iViscousFlux, iEulerXFlux
      public  iViscousFlux0D, iViscousFlux2D, iViscousFlux3D
      public  iEulerFlux0D, iEulerFlux3D

     interface iEulerFlux
         module procedure iEulerFlux0D, iEulerFlux3D
     end interface iEulerFlux

     interface iViscousFlux
         module procedure iViscousFlux0D, iViscousFlux2D, iViscousFlux3D
     end interface iViscousFlux
!
!     ========
      CONTAINS 
!     ========
!
!//////////////////////////////////////////////////////////////////////////////
!
!           INVISCID FLUXES
!           ---------------   
!
!//////////////////////////////////////////////////////////////////////////////
!
      pure subroutine iEulerFlux0D(Q, F)
         implicit none
         real(kind=RP), intent(in)   :: Q(1:NINC)
         real(kind=RP), intent(out)  :: F(1:NINC, 1:NDIM)
!
!        ---------------
!        Local variables
!        ---------------
!
         real(kind=RP) :: invRho

         invRho = 1.0_RP / Q(INSRHO)
!
!        X-Flux
!        ------         
         F(INSRHO , IX) = Q(INSRHOU)
         F(INSRHOU, IX) = invRho*Q(INSRHOU)*Q(INSRHOU) + Q(INSP)
         F(INSRHOV, IX) = invRho*Q(INSRHOU)*Q(INSRHOV)
         F(INSRHOW, IX) = invRho*Q(INSRHOU)*Q(INSRHOW)
         F(INSP   , IX) = thermodynamics % rho0c02 * invRho * Q(INSRHOU)
!
!        Y-Flux
!        ------
         F(INSRHO , IY) = Q(INSRHOV)
         F(INSRHOU, IY) = invRho*Q(INSRHOV)*Q(INSRHOU)
         F(INSRHOV, IY) = invRho*Q(INSRHOV)*Q(INSRHOV) + Q(INSP)
         F(INSRHOW, IY) = invRho*Q(INSRHOV)*Q(INSRHOW)
         F(INSP   , IY) = thermodynamics % rho0c02 * invRho * Q(INSRHOV)
!
!        Z-Flux
!        ------
         F(INSRHO ,IZ) = Q(INSRHOW)
         F(INSRHOU,IZ) = invRho*Q(INSRHOW)*Q(INSRHOU)
         F(INSRHOV,IZ) = invRho*Q(INSRHOW)*Q(INSRHOV)
         F(INSRHOW,IZ) = invRho*Q(INSRHOW)*Q(INSRHOW) + Q(INSP)
         F(INSP   ,IZ) = thermodynamics % rho0c02 * invRho * Q(INSRHOW)
      
      end subroutine iEulerFlux0D

      pure subroutine iEulerXFlux(Q, F)
         implicit none
         real(kind=RP), intent(in)   :: Q(1:NINC)
         real(kind=RP), intent(out)  :: F(1:NINC)
!
!        ---------------
!        Local variables
!        ---------------
!
         real(kind=RP) :: invRho

         invRho = 1.0_RP / Q(INSRHO)
!
!        X-Flux
!        ------         
         F(INSRHO ) = Q(INSRHOU)
         F(INSRHOU) = invRho*Q(INSRHOU)*Q(INSRHOU) + Q(INSP)
         F(INSRHOV) = invRho*Q(INSRHOU)*Q(INSRHOV)
         F(INSRHOW) = invRho*Q(INSRHOU)*Q(INSRHOW)
         F(INSP   ) = thermodynamics % rho0c02 * invRho * Q(INSRHOU)
      
      end subroutine iEulerXFlux

      pure subroutine iEulerFlux3D(N, Q, F)
         implicit none
         integer,       intent(in)  :: N(3)
         real(kind=RP), intent(in)  :: Q(1:NINC,0:N(1),0:N(2),0:N(3))
         real(kind=RP), intent(out) :: F(1:NINC,0:N(1),0:N(2),0:N(3),1:NDIM)
!
!        ---------------
!        Local variables
!        ---------------
!
         integer                 :: i, j, k
         real(kind=RP)           :: invRho

         do k = 0, N(3) ; do j = 0, N(2) ; do i = 0, N(1)
            invRho = 1.0_RP / Q(INSRHO,i,j,k)
!   
!           X-Flux
!           ------         
            F(INSRHO ,i,j,k, IX) = Q(INSRHOU,i,j,k)
            F(INSRHOU,i,j,k, IX) = invRho*Q(INSRHOU,i,j,k)*Q(INSRHOU,i,j,k) + Q(INSP,i,j,k)
            F(INSRHOV,i,j,k, IX) = invRho*Q(INSRHOU,i,j,k)*Q(INSRHOV,i,j,k)
            F(INSRHOW,i,j,k, IX) = invRho*Q(INSRHOU,i,j,k)*Q(INSRHOW,i,j,k)
            F(INSP   ,i,j,k, IX) = thermodynamics % rho0c02 * invRho * Q(INSRHOU,i,j,k)
!   
!           Y-Flux
!           ------
            F(INSRHO ,i,j,k, IY) = Q(INSRHOV,i,j,k)
            F(INSRHOU,i,j,k, IY) = invRho*Q(INSRHOV,i,j,k)*Q(INSRHOU,i,j,k)
            F(INSRHOV,i,j,k, IY) = invRho*Q(INSRHOV,i,j,k)*Q(INSRHOV,i,j,k) + Q(INSP,i,j,k)
            F(INSRHOW,i,j,k, IY) = invRho*Q(INSRHOV,i,j,k)*Q(INSRHOW,i,j,k)
            F(INSP   ,i,j,k, IY) = thermodynamics % rho0c02 * invRho * Q(INSRHOV,i,j,k)
!   
!           Z-Flux
!           ------
            F(INSRHO ,i,j,k, IZ) = Q(INSRHOW,i,j,k)
            F(INSRHOU,i,j,k, IZ) = invRho*Q(INSRHOW,i,j,k)*Q(INSRHOU,i,j,k)
            F(INSRHOV,i,j,k, IZ) = invRho*Q(INSRHOW,i,j,k)*Q(INSRHOV,i,j,k)
            F(INSRHOW,i,j,k, IZ) = invRho*Q(INSRHOW,i,j,k)*Q(INSRHOW,i,j,k) + Q(INSP,i,j,k)
            F(INSP   ,i,j,k, IZ) = thermodynamics % rho0c02 * invRho * Q(INSRHOW,i,j,k)
   
         end do   ; end do          ; end do

      end subroutine iEulerFlux3D
!
!//////////////////////////////////////////////////////////////////////////////////////////
!
!>        VISCOUS FLUXES
!         --------------
!
!//////////////////////////////////////////////////////////////////////////////////////////
!
      pure subroutine iViscousFlux0D(nEqn, nGradEqn, Q, U_x, U_y, U_z, mu, beta, kappa, F)
         implicit none
         integer,       intent(in)  :: nEqn
         integer,       intent(in)  :: nGradEqn
         real(kind=RP), intent(in)  :: Q   (1:nEqn     )
         real(kind=RP), intent(in)  :: U_x (1:nGradEqn)
         real(kind=RP), intent(in)  :: U_y (1:nGradEqn)
         real(kind=RP), intent(in)  :: U_z (1:nGradEqn)
         real(kind=RP), intent(in)  :: mu
         real(kind=RP), intent(in)  :: beta
         real(kind=RP), intent(in)  :: kappa
         real(kind=RP), intent(out) :: F(1:nEqn, 1:NDIM)
!
!        ---------------
!        Local variables
!        ---------------
!
         real(kind=RP) :: gradV_x(NDIM), gradV_y(NDIM), gradV_z(NDIM), invRho, uDivRho(NDIM)

         invRho = 1.0_RP / Q(INSRHO)

         uDivRho = Q(INSRHOU:INSRHOW) * invRho * invRho

         gradV_x = invRho * U_x(INSRHOU:INSRHOW) - uDivRho * U_x(INSRHO)
         gradV_y = invRho * U_y(INSRHOU:INSRHOW) - uDivRho * U_y(INSRHO)
         gradV_z = invRho * U_z(INSRHOU:INSRHOW) - uDivRho * U_z(INSRHO)
         

         F(INSRHO,IX)  = 0.0_RP
         F(INSRHOU,IX) = mu * gradV_x(IX)
         F(INSRHOV,IX) = mu * gradV_x(IY)
         F(INSRHOW,IX) = mu * gradV_x(IZ)
         F(INSP,IX) = 0.0_RP

         F(INSRHO,IY)  = 0.0_RP
         F(INSRHOU,IY) = mu * gradV_y(IX)
         F(INSRHOV,IY) = mu * gradV_y(IY)
         F(INSRHOW,IY) = mu * gradV_y(IZ)
         F(INSP,IY) = 0.0_RP

         F(INSRHO,IZ)  = 0.0_RP
         F(INSRHOU,IZ) = mu * gradV_z(IX)
         F(INSRHOV,IZ) = mu * gradV_z(IY)
         F(INSRHOW,IZ) = mu * gradV_z(IZ)
         F(INSP,IZ) = 0.0_RP

      end subroutine iViscousFlux0D

      pure subroutine iViscousFlux2D( nEqn, nGradEqn, N, Q, U_x, U_y, U_z, mu, beta, kappa, F)
         implicit none
         integer,       intent(in)  :: nEqn
         integer,       intent(in)  :: nGradEqn
         integer         , intent(in)  :: N(2)
         real(kind=RP),    intent(in)  :: Q  (1:nEqn, 0:N(1), 0:N(2))
         real(kind=RP),    intent(in)  :: U_x(1:nGradEqn, 0:N(1), 0:N(2) )
         real(kind=RP),    intent(in)  :: U_y(1:nGradEqn, 0:N(1), 0:N(2) )
         real(kind=RP),    intent(in)  :: U_z(1:nGradEqn, 0:N(1), 0:N(2) )
         real(kind=RP),    intent(in)  :: mu  (0:N(1), 0:N(2))
         real(kind=RP),    intent(in)  :: beta(0:N(1), 0:N(2))
         real(kind=RP),    intent(in)  :: kappa(0:N(1), 0:N(2))
         real(kind=RP),    intent(out) :: F   (1:nEqn, 1:NDIM, 0:N(1), 0:N(2))
!
!        ---------------
!        Local variables
!        ---------------
!
         integer       :: i , j
         real(kind=RP) :: gradV_x(NDIM), gradV_y(NDIM), gradV_z(NDIM), invRho, uDivRho(NDIM)

         do j = 0, N(2) ; do i = 0, N(1)
            invRho = 1.0_RP / Q(INSRHO,i,j)
   
            uDivRho = Q(INSRHOU:INSRHOW,i,j) * invRho * invRho
   
            gradV_x = invRho * U_x(INSRHOU:INSRHOW,i,j) - uDivRho * U_x(INSRHO,i,j)
            gradV_y = invRho * U_y(INSRHOU:INSRHOW,i,j) - uDivRho * U_y(INSRHO,i,j)
            gradV_z = invRho * U_z(INSRHOU:INSRHOW,i,j) - uDivRho * U_z(INSRHO,i,j)
            
   
            F(INSRHO,IX,i,j)  = 0.0_RP
            F(INSRHOU,IX,i,j) = mu(i,j) * gradV_x(IX)
            F(INSRHOV,IX,i,j) = mu(i,j) * gradV_x(IY)
            F(INSRHOW,IX,i,j) = mu(i,j) * gradV_x(IZ)
            F(INSP,IX,i,j) = 0.0_RP
   
            F(INSRHO,IY,i,j)  = 0.0_RP
            F(INSRHOU,IY,i,j) = mu(i,j) * gradV_y(IX)
            F(INSRHOV,IY,i,j) = mu(i,j) * gradV_y(IY)
            F(INSRHOW,IY,i,j) = mu(i,j) * gradV_y(IZ)
            F(INSP,IY,i,j) = 0.0_RP
   
            F(INSRHO,IZ,i,j)  = 0.0_RP
            F(INSRHOU,IZ,i,j) = mu(i,j) * gradV_z(IX)
            F(INSRHOV,IZ,i,j) = mu(i,j) * gradV_z(IY)
            F(INSRHOW,IZ,i,j) = mu(i,j) * gradV_z(IZ)
            F(INSP,IZ,i,j) = 0.0_RP
   
         end do    ; end do

      end subroutine iViscousFlux2D

      pure subroutine iViscousFlux3D( nEqn, nGradEqn, N, Q, U_x, U_y, U_z, mu, beta, kappa, F)
         implicit none
         integer,       intent(in)  :: nEqn
         integer,       intent(in)  :: nGradEqn
         integer         , intent(in)  :: N(3)
         real(kind=RP),    intent(in)  :: Q  (1:nEqn, 0:N(1), 0:N(2), 0:N(3))
         real(kind=RP),    intent(in)  :: U_x(1:nGradEqn, 0:N(1), 0:N(2), 0:N(3) )
         real(kind=RP),    intent(in)  :: U_y(1:nGradEqn, 0:N(1), 0:N(2), 0:N(3) )
         real(kind=RP),    intent(in)  :: U_z(1:nGradEqn, 0:N(1), 0:N(2), 0:N(3) )
         real(kind=RP),    intent(in)  :: mu  (0:N(1), 0:N(2), 0:N(3))
         real(kind=RP),    intent(in)  :: beta(0:N(1), 0:N(2), 0:N(3))
         real(kind=RP),    intent(in)  :: kappa(0:N(1), 0:N(2), 0:N(3))
         real(kind=RP),    intent(out) :: F   (1:nEqn, 0:N(1), 0:N(2), 0:N(3), 1:NDIM )
!
!        ---------------
!        Local variables
!        ---------------
!
         integer       :: i , j , k
         real(kind=RP) :: gradV_x(NDIM), gradV_y(NDIM), gradV_z(NDIM), invRho, uDivRho(NDIM)

         do k = 0, N(3) ; do j = 0, N(2) ; do i = 0, N(1)
   
            invRho = 1.0_RP / Q(INSRHO,i,j,k)
   
            uDivRho = Q(INSRHOU:INSRHOW,i,j,k) * invRho * invRho
   
            gradV_x = invRho * U_x(INSRHOU:INSRHOW,i,j,k) - uDivRho * U_x(INSRHO,i,j,k)
            gradV_y = invRho * U_y(INSRHOU:INSRHOW,i,j,k) - uDivRho * U_y(INSRHO,i,j,k)
            gradV_z = invRho * U_z(INSRHOU:INSRHOW,i,j,k) - uDivRho * U_z(INSRHO,i,j,k)
            
   
            F(INSRHO,i,j,k,IX)  = 0.0_RP
            F(INSRHOU,i,j,k,IX) = mu(i,j,k) * gradV_x(IX)
            F(INSRHOV,i,j,k,IX) = mu(i,j,k) * gradV_x(IY)
            F(INSRHOW,i,j,k,IX) = mu(i,j,k) * gradV_x(IZ)
            F(INSP,i,j,k,IX) = 0.0_RP
   
            F(INSRHO,i,j,k,IY)  = 0.0_RP
            F(INSRHOU,i,j,k,IY) = mu(i,j,k) * gradV_y(IX)
            F(INSRHOV,i,j,k,IY) = mu(i,j,k) * gradV_y(IY)
            F(INSRHOW,i,j,k,IY) = mu(i,j,k) * gradV_y(IZ)
            F(INSP,i,j,k,IY) = 0.0_RP
   
            F(INSRHO,i,j,k,IZ)  = 0.0_RP
            F(INSRHOU,i,j,k,IZ) = mu(i,j,k) * gradV_z(IX)
            F(INSRHOV,i,j,k,IZ) = mu(i,j,k) * gradV_z(IY)
            F(INSRHOW,i,j,k,IZ) = mu(i,j,k) * gradV_z(IZ)
            F(INSP,i,j,k,IZ) = 0.0_RP
   
         end do      ; end do    ; end do

      end subroutine iViscousFlux3D

   END Module Physics_iNS
!@mark -
!
! /////////////////////////////////////////////////////////////////////
!
!----------------------------------------------------------------------
!! This routine returns the maximum eigenvalues for the Euler equations 
!! for the given solution value in each spatial direction. 
!! These are to be used to compute the local time step.
!----------------------------------------------------------------------
!
      SUBROUTINE ComputeEigenvaluesForState( Q, eigen )
      
      USE SMConstants
      USE PhysicsStorage_iNS
      use FluidData_iNS,          only: Thermodynamics
      IMPLICIT NONE
!
!     ---------
!     Arguments
!     ---------
!
      REAL(KIND=Rp), DIMENSION(NINC) :: Q
      REAL(KIND=Rp), DIMENSION(3)     :: eigen
!
!     ---------------
!     Local Variables
!     ---------------
!
      REAL(KIND=Rp) :: u, v, w, p, a
!      
      u = ABS( Q(INSRHOU)/Q(INSRHO) )
      v = ABS( Q(INSRHOV)/Q(INSRHO) )
      w = ABS( Q(INSRHOW)/Q(INSRHO) )
      a = sqrt(max(max(u,v),w)**2 + 4.0_RP * thermodynamics % rho0c02/Q(INSRHO))
      
      eigen(1) = u + a
      eigen(2) = v + a
      eigen(3) = w + a
      
      END SUBROUTINE ComputeEigenvaluesForState
