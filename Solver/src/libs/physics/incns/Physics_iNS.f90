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
      pure subroutine iEulerFlux0D(Q, F, rho_)
         implicit none
         real(kind=RP), intent(in)   :: Q(1:NCONS)
         real(kind=RP), intent(out)  :: F(1:NCONS, 1:NDIM)
         real(kind=RP), intent(in), optional :: rho_
!
!        ---------------
!        Local variables
!        ---------------
!
         real(kind=RP) :: rho, invRho

         rho = (Q(INSRHO))

         invRho = 1.0_RP / rho
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
         real(kind=RP), intent(in)   :: Q(1:NCONS)
         real(kind=RP), intent(out)  :: F(1:NCONS)
!
!        ---------------
!        Local variables
!        ---------------
!
         real(kind=RP) :: rho, invRho

         rho = (Q(INSRHO))
      
         invRho = 1.0_RP / rho
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
         real(kind=RP), intent(in)  :: Q(1:NCONS,0:N(1),0:N(2),0:N(3))
         real(kind=RP), intent(out) :: F(1:NCONS,0:N(1),0:N(2),0:N(3),1:NDIM)
!
!        ---------------
!        Local variables
!        ---------------
!
         integer                 :: i, j, k
         real(kind=RP)           :: rho, invRho

         do k = 0, N(3) ; do j = 0, N(2) ; do i = 0, N(1)
            rho = (Q(INSRHO,i,j,k))
            invRho = 1.0_RP / rho
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

         F(INSRHO,IX)  = 0.0_RP
         F(INSRHOU,IX) = 2.0_RP * mu * U_x(INSRHOU)
         F(INSRHOV,IX) = mu * (U_x(INSRHOV) + U_y(INSRHOU))
         F(INSRHOW,IX) = mu * (U_x(INSRHOW) + U_z(INSRHOU))
         F(INSP,IX) = 0.0_RP

         F(INSRHO,IY)  = 0.0_RP
         F(INSRHOU,IY) = F(INSRHOV,IX)
         F(INSRHOV,IY) = 2.0_RP * mu * U_y(INSRHOV)
         F(INSRHOW,IY) = mu * (U_y(INSRHOW) + U_z(INSRHOV))
         F(INSP,IY) = 0.0_RP

         F(INSRHO,IZ)  = 0.0_RP
         F(INSRHOU,IZ) = F(INSRHOW,IX)
         F(INSRHOV,IZ) = F(INSRHOW,IY)
         F(INSRHOW,IZ) = 2.0_RP * mu * U_z(INSRHOW)
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

         do j = 0, N(2) ; do i = 0, N(1)
            F(INSRHO,IX,i,j)  = 0.0_RP
            F(INSRHOU,IX,i,j) = 2.0_RP * mu(i,j) * U_x(INSRHOU,i,j)
            F(INSRHOV,IX,i,j) = mu(i,j) * (U_x(INSRHOV,i,j) + U_y(INSRHOU,i,j))
            F(INSRHOW,IX,i,j) = mu(i,j) * (U_x(INSRHOW,i,j) + U_z(INSRHOU,i,j))
            F(INSP,IX,i,j) = 0.0_RP
   
            F(INSRHO,IY,i,j)  = 0.0_RP
            F(INSRHOU,IY,i,j) = F(INSRHOV,IX,i,j)
            F(INSRHOV,IY,i,j) = 2.0_RP * mu(i,j) * U_y(INSRHOV,i,j)
            F(INSRHOW,IY,i,j) = mu(i,j) * (U_y(INSRHOW,i,j) + U_z(INSRHOV,i,j))
            F(INSP,IY,i,j) = 0.0_RP
   
            F(INSRHO,IZ,i,j)  = 0.0_RP
            F(INSRHOU,IZ,i,j) = F(INSRHOW,IX,i,j)
            F(INSRHOV,IZ,i,j) = F(INSRHOW,IY,i,j)
            F(INSRHOW,IZ,i,j) = 2.0_RP * mu(i,j) * U_z(INSRHOW,i,j)
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

         do k = 0, N(3) ; do j = 0, N(2) ; do i = 0, N(1)
            F(INSRHO,i,j,k,IX)  = 0.0_RP
            F(INSRHOU,i,j,k,IX) = 2.0_RP * mu(i,j,k) * U_x(INSRHOU,i,j,k)
            F(INSRHOV,i,j,k,IX) = mu(i,j,k) * (U_x(INSRHOV,i,j,k) + U_y(INSRHOU,i,j,k))
            F(INSRHOW,i,j,k,IX) = mu(i,j,k) * (U_x(INSRHOW,i,j,k) + U_z(INSRHOU,i,j,k))
            F(INSP,i,j,k,IX) = 0.0_RP
   
            F(INSRHO,i,j,k,IY)  = 0.0_RP
            F(INSRHOU,i,j,k,IY) = F(INSRHOV,i,j,k,IX)
            F(INSRHOV,i,j,k,IY) = 2.0_RP * mu(i,j,k) * U_y(INSRHOV,i,j,k)
            F(INSRHOW,i,j,k,IY) = mu(i,j,k) * (U_y(INSRHOW,i,j,k) + U_z(INSRHOV,i,j,k))
            F(INSP,i,j,k,IY) = 0.0_RP
   
            F(INSRHO,i,j,k,IZ)  = 0.0_RP
            F(INSRHOU,i,j,k,IZ) = F(INSRHOW,i,j,k,IX)
            F(INSRHOV,i,j,k,IZ) = F(INSRHOW,i,j,k,IY)
            F(INSRHOW,i,j,k,IZ) = 2.0_RP * mu(i,j,k) * U_z(INSRHOW,i,j,k)
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
      REAL(KIND=Rp), DIMENSION(NCONS) :: Q
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
      a = sqrt(u*u+v*v+w*w + 4.0_RP * thermodynamics % rho0c02/Q(INSRHO))
      
      eigen(1) = 0.5_RP * (u + a)
      eigen(2) = 0.5_RP * (v + a)
      eigen(3) = 0.5_RP * (w + a)
      
      END SUBROUTINE ComputeEigenvaluesForState
