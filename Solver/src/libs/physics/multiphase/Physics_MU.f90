!
!//////////////////////////////////////////////////////
!
!   @File:    Physics_MU.f90
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
   module Physics_MU
!  **************
!
      use SMConstants
      use PhysicsStorage_MU
      use VariableConversion_MU
      use FluidData_MU
      implicit none

      private
      public  mEulerFlux, mViscousFlux, mEulerXFlux
      public  mViscousFlux0D, mViscousFlux2D, mViscousFlux3D
      public  mEulerFlux0D, mEulerFlux3D

     interface mEulerFlux
         module procedure mEulerFlux0D, mEulerFlux3D
     end interface mEulerFlux

     interface mViscousFlux
         module procedure mViscousFlux0D, mViscousFlux2D, mViscousFlux3D
     end interface mViscousFlux
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
      pure subroutine mEulerFlux0D(Q, F, rho_)
         implicit none
         real(kind=RP), intent(in)   :: Q(1:NCONS)
         real(kind=RP), intent(out)  :: F(1:NCONS, 1:NDIM)
         real(kind=RP), intent(in), optional :: rho_
!
!        ---------------
!        Local variables
!        ---------------
!
         real(kind=RP) :: invSqrtRho

         invSqrtRho = 1.0_RP / sqrt(rho_)
!
!        X-Flux
!        ------         
         F(IMC     , IX) = Q(IMC)*invSqrtRho*Q(IMSQRHOU)
         F(IMSQRHOU, IX) = 0.5_RP*Q(IMSQRHOU)*Q(IMSQRHOU) + Q(IMP)
         F(IMSQRHOV, IX) = 0.5_RP*Q(IMSQRHOU)*Q(IMSQRHOV)
         F(IMSQRHOW, IX) = 0.5_RP*Q(IMSQRHOU)*Q(IMSQRHOW)
         F(IMP     , IX) = thermodynamics % rho0c02 * invSqrtRho * Q(IMSQRHOU)
!
!        Y-Flux
!        ------
         F(IMC     , IY) = Q(IMC)*invSqrtRho*Q(IMSQRHOV)
         F(IMSQRHOU, IY) = 0.5_RP*Q(IMSQRHOV)*Q(IMSQRHOU)
         F(IMSQRHOV, IY) = 0.5_RP*Q(IMSQRHOV)*Q(IMSQRHOV) + Q(IMP)
         F(IMSQRHOW, IY) = 0.5_RP*Q(IMSQRHOV)*Q(IMSQRHOW)
         F(IMP     , IY) = thermodynamics % rho0c02 * invSqrtRho * Q(IMSQRHOV)
!
!        Z-Flux
!        ------
         F(IMC     , IZ) = Q(IMC)*invSqrtRho*Q(IMSQRHOW)
         F(IMSQRHOU, IZ) = 0.5_RP*Q(IMSQRHOW)*Q(IMSQRHOU)
         F(IMSQRHOV, IZ) = 0.5_RP*Q(IMSQRHOW)*Q(IMSQRHOV) 
         F(IMSQRHOW, IZ) = 0.5_RP*Q(IMSQRHOW)*Q(IMSQRHOW) + Q(IMP)
         F(IMP     , IZ) = thermodynamics % rho0c02 * invSqrtRho * Q(IMSQRHOW)
      
      end subroutine mEulerFlux0D

      pure subroutine mEulerXFlux(Q, F, rho)
         implicit none
         real(kind=RP), intent(in)   :: Q(1:NCONS)
         real(kind=RP), intent(out)  :: F(1:NCONS)
         real(kind=RP), intent(in)   :: rho
!
!        ---------------
!        Local variables
!        ---------------
!
         real(kind=RP) :: invSqrtRho

         invSqrtRho = 1.0_RP / sqrt(rho)
!
!        X-Flux
!        ------         
         F(IMC)      = Q(IMC)*invSqrtRho*Q(IMSQRHOU)
         F(IMSQRHOU) = 0.5_RP*Q(IMSQRHOU)*Q(IMSQRHOU) + Q(IMP)
         F(IMSQRHOV) = 0.5_RP*Q(IMSQRHOU)*Q(IMSQRHOV)
         F(IMSQRHOW) = 0.5_RP*Q(IMSQRHOU)*Q(IMSQRHOW)
         F(IMP)      = thermodynamics % rho0c02 * invSqrtRho * Q(IMSQRHOU)
      
      end subroutine mEulerXFlux

      pure subroutine mEulerFlux3D(N, Q, F, rho_)
         implicit none
         integer,       intent(in)  :: N(3)
         real(kind=RP), intent(in)  :: Q(1:NCONS,0:N(1),0:N(2),0:N(3))
         real(kind=RP), intent(out) :: F(1:NCONS,0:N(1),0:N(2),0:N(3),1:NDIM)
         real(kind=RP), intent(in), optional :: rho_(0:N(1), 0:N(2), 0:N(3))
!
!        ---------------
!        Local variables
!        ---------------
!
         integer                 :: i, j, k
         real(kind=RP)           :: invSqrtRho

         do k = 0, N(3) ; do j = 0, N(2) ; do i = 0, N(1)
            invSqrtRho = 1.0_RP / sqrt(rho_(i,j,k))
!   
!           X-Flux
!           ------         
            F(IMC     ,i,j,k, IX) = Q(IMC,i,j,k)*invSqrtRho*Q(IMSQRHOU,i,j,k)
            F(IMSQRHOU,i,j,k, IX) = 0.5_RP*Q(IMSQRHOU,i,j,k)*Q(IMSQRHOU,i,j,k) + Q(IMP,i,j,k)
            F(IMSQRHOV,i,j,k, IX) = 0.5_RP*Q(IMSQRHOU,i,j,k)*Q(IMSQRHOV,i,j,k)
            F(IMSQRHOW,i,j,k, IX) = 0.5_RP*Q(IMSQRHOU,i,j,k)*Q(IMSQRHOW,i,j,k)
            F(IMP     ,i,j,k, IX) = thermodynamics % rho0c02 * invSqrtRho * Q(IMSQRHOU,i,j,k)
!   
!           Y-Flux
!           ------
            F(IMC     ,i,j,k, IY) = Q(IMC,i,j,k)*invSqrtRho*Q(IMSQRHOV,i,j,k)
            F(IMSQRHOU,i,j,k, IY) = 0.5_RP*Q(IMSQRHOV,i,j,k)*Q(IMSQRHOU,i,j,k)
            F(IMSQRHOV,i,j,k, IY) = 0.5_RP*Q(IMSQRHOV,i,j,k)*Q(IMSQRHOV,i,j,k) + Q(IMP,i,j,k)
            F(IMSQRHOW,i,j,k, IY) = 0.5_RP*Q(IMSQRHOV,i,j,k)*Q(IMSQRHOW,i,j,k)
            F(IMP     ,i,j,k, IY) = thermodynamics % rho0c02 * invSqrtRho * Q(IMSQRHOV,i,j,k)
!   
!           Z-Flux
!           ------
            F(IMC     ,i,j,k, IZ) = Q(IMC,i,j,k)*invSqrtRho*Q(IMSQRHOW,i,j,k)
            F(IMSQRHOU,i,j,k, IZ) = 0.5_RP*Q(IMSQRHOW,i,j,k)*Q(IMSQRHOU,i,j,k)
            F(IMSQRHOV,i,j,k, IZ) = 0.5_RP*Q(IMSQRHOW,i,j,k)*Q(IMSQRHOV,i,j,k) 
            F(IMSQRHOW,i,j,k, IZ) = 0.5_RP*Q(IMSQRHOW,i,j,k)*Q(IMSQRHOW,i,j,k) + Q(IMP,i,j,k)
            F(IMP     ,i,j,k, IZ) = thermodynamics % rho0c02 * invSqrtRho * Q(IMSQRHOW,i,j,k)

         end do   ; end do          ; end do

      end subroutine mEulerFlux3D
!
!//////////////////////////////////////////////////////////////////////////////////////////
!
!>        VISCOUS FLUXES
!         --------------
!
!//////////////////////////////////////////////////////////////////////////////////////////
!
      pure subroutine mViscousFlux0D(nEqn, nGradEqn, Q, U_x, U_y, U_z, mu, beta, kappa, F)
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

         F(IGMU,IX)  = 0.0_RP
         F(IGU,IX) = 2.0_RP * mu * U_x(IGU)
         F(IGV,IX) = mu * (U_x(IGV) + U_y(IGU))
         F(IGW,IX) = mu * (U_x(IGW) + U_z(IGU))
         F(IGP,IX) = 0.0_RP

         F(IGMU,IY)  = 0.0_RP
         F(IGU,IY) = F(IGV,IX)
         F(IGV,IY) = 2.0_RP * mu * U_y(IGV)
         F(IGW,IY) = mu * (U_y(IGW) + U_z(IGV))
         F(IGP,IY) = 0.0_RP

         F(IGMU,IZ)  = 0.0_RP
         F(IGU,IZ) = F(IGW,IX)
         F(IGV,IZ) = F(IGW,IY)
         F(IGW,IZ) = 2.0_RP * mu * U_z(IGW)
         F(IGP,IZ) = 0.0_RP

      end subroutine mViscousFlux0D

      pure subroutine mViscousFlux2D( nEqn, nGradEqn, N, Q, U_x, U_y, U_z, mu, beta, kappa, F)
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
            F(IGMU,IX,i,j)  = 0.0_RP
            F(IGU,IX,i,j) = 2.0_RP * mu(i,j) * U_x(IGU,i,j)
            F(IGV,IX,i,j) = mu(i,j) * (U_x(IGV,i,j) + U_y(IGU,i,j))
            F(IGW,IX,i,j) = mu(i,j) * (U_x(IGW,i,j) + U_z(IGU,i,j))
            F(IGP,IX,i,j) = 0.0_RP
   
            F(IGMU,IY,i,j)  = 0.0_RP
            F(IGU,IY,i,j) = F(IGV,IX,i,j)
            F(IGV,IY,i,j) = 2.0_RP * mu(i,j) * U_y(IGV,i,j)
            F(IGW,IY,i,j) = mu(i,j) * (U_y(IGW,i,j) + U_z(IGV,i,j))
            F(IGP,IY,i,j) = 0.0_RP
   
            F(IGMU,IZ,i,j)  = 0.0_RP
            F(IGU,IZ,i,j) = F(IGW,IX,i,j)
            F(IGV,IZ,i,j) = F(IGW,IY,i,j)
            F(IGW,IZ,i,j) = 2.0_RP * mu(i,j) * U_z(IGW,i,j)
            F(IGP,IZ,i,j) = 0.0_RP
         end do    ; end do

      end subroutine mViscousFlux2D

      pure subroutine mViscousFlux3D( nEqn, nGradEqn, N, Q, U_x, U_y, U_z, mu, beta, kappa, F)
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
            F(IGMU,i,j,k,IX)  = 0.0_RP
            F(IGU,i,j,k,IX) = 2.0_RP * mu(i,j,k) * U_x(IGU,i,j,k)
            F(IGV,i,j,k,IX) = mu(i,j,k) * (U_x(IGV,i,j,k) + U_y(IGU,i,j,k))
            F(IGW,i,j,k,IX) = mu(i,j,k) * (U_x(IGW,i,j,k) + U_z(IGU,i,j,k))
            F(IGP,i,j,k,IX) = 0.0_RP
   
            F(IGMU,i,j,k,IY)  = 0.0_RP
            F(IGU,i,j,k,IY) = F(IGV,i,j,k,IX)
            F(IGV,i,j,k,IY) = 2.0_RP * mu(i,j,k) * U_y(IGV,i,j,k)
            F(IGW,i,j,k,IY) = mu(i,j,k) * (U_y(IGW,i,j,k) + U_z(IGV,i,j,k))
            F(IGP,i,j,k,IY) = 0.0_RP
   
            F(IGMU,i,j,k,IZ)  = 0.0_RP
            F(IGU,i,j,k,IZ) = F(IGW,i,j,k,IX)
            F(IGV,i,j,k,IZ) = F(IGW,i,j,k,IY)
            F(IGW,i,j,k,IZ) = 2.0_RP * mu(i,j,k) * U_z(IGW,i,j,k)
            F(IGP,i,j,k,IZ) = 0.0_RP
         end do      ; end do    ; end do

      end subroutine mViscousFlux3D

   END Module Physics_MU
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
      USE PhysicsStorage_MU
      use FluidData_MU,          only: Thermodynamics
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
print*, "Get eigenvalues!!"
errorMessage(STD_OUT)
stop
!      
!      u = ABS( Q(IMSQRHOU)/Q(IMSQRHO) )
!      v = ABS( Q(IMSQRHOV)/Q(IMSQRHO) )
!      w = ABS( Q(IMSQRHOW)/Q(IMSQRHO) )
!      a = sqrt(u*u+v*v+w*w + 4.0_RP * thermodynamics % rho0c02/Q(IMSQRHO))
      
!      eigen(1) = 0.5_RP * (u + a)
!      eigen(2) = 0.5_RP * (v + a)
!      eigen(3) = 0.5_RP * (w + a)
       eigen = 0.0_RP
      
      END SUBROUTINE ComputeEigenvaluesForState
