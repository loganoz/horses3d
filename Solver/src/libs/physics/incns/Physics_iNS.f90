!
!//////////////////////////////////////////////////////
!
!   @File:    Physics_iNS.f90
!   @Author:  Juan Manzanero (juan.manzanero@upm.es)
!   @Created: Tue Jun 19 17:39:26 2018
!   @Last revision date: Sat Jun 23 10:20:35 2018
!   @Last revision author: Juan Manzanero (juan.manzanero@upm.es)
!   @Last revision commit: fce351220409e80ce5df1949249c2b870dd847aa
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
      public  iEulerFlux, iViscousFlux
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
         real(kind=RP)           :: u , v , w , p

!
!        X-Flux
!        ------         
         F(INSRHO , IX ) = Q(INSRHO)*Q(INSU) 
         F(INSU, IX ) = Q(INSRHO)*Q(INSU)*Q(INSU) + Q(INSP)
         F(INSV, IX ) = Q(INSRHO)*Q(INSU)*Q(INSV)
         F(INSW, IX ) = Q(INSRHO)*Q(INSU)*Q(INSW)
         F(INSP, IX ) = thermodynamics % rho0c02 * Q(INSU)
!
!        Y-Flux
!        ------
         F(INSRHO , IY ) = Q(INSRHO)*Q(INSV)
         F(INSU ,IY ) = Q(INSRHO)*Q(INSV)*Q(INSU)
         F(INSV ,IY ) = Q(INSRHO)*Q(INSV)*Q(INSV) + Q(INSP)
         F(INSW ,IY ) = Q(INSRHO)*Q(INSV)*Q(INSW)
         F(INSP ,IY ) = thermodynamics % rho0c02 * Q(INSV)
!
!        Z-Flux
!        ------
         F(INSRHO ,IZ) = Q(INSRHO)*Q(INSW)
         F(INSU,IZ) = Q(INSRHO)*Q(INSW)*Q(INSU)
         F(INSV,IZ) = Q(INSRHO)*Q(INSW)*Q(INSV)
         F(INSW,IZ) = Q(INSRHO)*Q(INSW)*Q(INSW) + Q(INSP)
         F(INSP,IZ) = thermodynamics % rho0c02 * Q(INSW)
      
      end subroutine iEulerFlux0D

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

         do k = 0, N(3) ; do j = 0, N(2) ; do i = 0, N(1)
!   
!           X-Flux
!           ------         
            F(INSRHO ,i,j,k, IX) = Q(INSRHO,i,j,k)*Q(INSU,i,j,k) 
            F(INSU,i,j,k, IX) = Q(INSRHO,i,j,k)*Q(INSU,i,j,k)*Q(INSU,i,j,k) + Q(INSP,i,j,k)
            F(INSV,i,j,k, IX) = Q(INSRHO,i,j,k)*Q(INSU,i,j,k)*Q(INSV,i,j,k)
            F(INSW,i,j,k, IX) = Q(INSRHO,i,j,k)*Q(INSU,i,j,k)*Q(INSW,i,j,k)
            F(INSP,i,j,k, IX) = thermodynamics % rho0c02 * Q(INSU,i,j,k)
!   
!           Y-Flux
!           ------
            F(INSRHO ,i,j,k, IY) = Q(INSRHO,i,j,k)*Q(INSV,i,j,k)
            F(INSU ,i,j,k,IY) = Q(INSRHO,i,j,k)*Q(INSV,i,j,k)*Q(INSU,i,j,k)
            F(INSV ,i,j,k,IY) = Q(INSRHO,i,j,k)*Q(INSV,i,j,k)*Q(INSV,i,j,k) + Q(INSP,i,j,k)
            F(INSW ,i,j,k,IY) = Q(INSRHO,i,j,k)*Q(INSV,i,j,k)*Q(INSW,i,j,k)
            F(INSP ,i,j,k,IY) = thermodynamics % rho0c02 * Q(INSV,i,j,k)
!   
!           Z-Flux
!           ------
            F(INSRHO ,i,j,k,IZ) = Q(INSRHO,i,j,k)*Q(INSW,i,j,k)
            F(INSU,i,j,k,IZ) = Q(INSRHO,i,j,k)*Q(INSW,i,j,k)*Q(INSU,i,j,k)
            F(INSV,i,j,k,IZ) = Q(INSRHO,i,j,k)*Q(INSW,i,j,k)*Q(INSV,i,j,k)
            F(INSW,i,j,k,IZ) = Q(INSRHO,i,j,k)*Q(INSW,i,j,k)*Q(INSW,i,j,k) + Q(INSP,i,j,k)
            F(INSP,i,j,k,IZ) = thermodynamics % rho0c02 * Q(INSW,i,j,k)

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
      pure subroutine iViscousFlux0D(nEqn, nGradEqn, Q, U_x, U_y, U_z, mu, kappa, F)
         implicit none
         integer,       intent(in)  :: nEqn
         integer,       intent(in)  :: nGradEqn
         real(kind=RP), intent(in)  :: Q   (1:nEqn     )
         real(kind=RP), intent(in)  :: U_x (1:nGradEqn)
         real(kind=RP), intent(in)  :: U_y (1:nGradEqn)
         real(kind=RP), intent(in)  :: U_z (1:nGradEqn)
         real(kind=RP), intent(in)  :: mu
         real(kind=RP), intent(in)  :: kappa
         real(kind=RP), intent(out) :: F(1:nEqn, 1:NDIM)

         F(INSRHO,IX)  = 0.0_RP
         F(INSU,IX) = mu * U_x(INSU)
         F(INSV,IX) = mu * U_x(INSV)
         F(INSW,IX) = mu * U_x(INSW)
         F(INSP,IX) = 0.0_RP

         F(INSRHO,IY)  = 0.0_RP
         F(INSU,IY) = mu * U_y(INSU)
         F(INSV,IY) = mu * U_y(INSV)
         F(INSW,IY) = mu * U_y(INSW)
         F(INSP,IY) = 0.0_RP

         F(INSRHO,IZ)  = 0.0_RP
         F(INSU,IZ) = mu * U_z(INSU)
         F(INSV,IZ) = mu * U_z(INSV)
         F(INSW,IZ) = mu * U_z(INSW)
         F(INSP,IZ) = 0.0_RP

      end subroutine iViscousFlux0D

      pure subroutine iViscousFlux2D( nEqn, nGradEqn, N, Q, U_x, U_y, U_z, mu, kappa, F)
         implicit none
         integer,       intent(in)  :: nEqn
         integer,       intent(in)  :: nGradEqn
         integer         , intent(in)  :: N(2)
         real(kind=RP),    intent(in)  :: Q  (1:nEqn, 0:N(1), 0:N(2))
         real(kind=RP),    intent(in)  :: U_x(1:nGradEqn, 0:N(1), 0:N(2) )
         real(kind=RP),    intent(in)  :: U_y(1:nGradEqn, 0:N(1), 0:N(2) )
         real(kind=RP),    intent(in)  :: U_z(1:nGradEqn, 0:N(1), 0:N(2) )
         real(kind=RP),    intent(in)  :: mu  (0:N(1), 0:N(2))
         real(kind=RP),    intent(in)  :: kappa(0:N(1), 0:N(2))
         real(kind=RP),    intent(out) :: F   (1:nEqn, 1:NDIM, 0:N(1), 0:N(2))
!
!        ---------------
!        Local variables
!        ---------------
!
         integer       :: i , j

         do j = 0, N(2) ; do i = 0, N(1)
            F(INSRHO ,IX,i,j) = 0.0_RP
            F(INSU,IX,i,j)    = mu(i,j) * U_x(INSU,i,j)
            F(INSV,IX,i,j)    = mu(i,j) * U_x(INSV,i,j)
            F(INSW,IX,i,j)    = mu(i,j) * U_x(INSW,i,j) 
            F(INSP,IX,i,j)    = 0.0_RP
   
            F(INSRHO ,IY,i,j) = 0.0_RP
            F(INSU,IY,i,j)    = mu(i,j) * U_y(INSU,i,j)
            F(INSV,IY,i,j)    = mu(i,j) * U_y(INSV,i,j)
            F(INSW,IY,i,j)    = mu(i,j) * U_y(INSW,i,j) 
            F(INSP,IY,i,j)    = 0.0_RP
   
            F(INSRHO ,IZ,i,j) = 0.0_RP
            F(INSU,IZ,i,j)    = mu(i,j) * U_z(INSU,i,j)
            F(INSV,IZ,i,j)    = mu(i,j) * U_z(INSV,i,j)
            F(INSW,IZ,i,j)    = mu(i,j) * U_z(INSW,i,j) 
            F(INSP,IZ,i,j)    = 0.0_RP
   
         end do    ; end do

      end subroutine iViscousFlux2D

      pure subroutine iViscousFlux3D( nEqn, nGradEqn, N, Q, U_x, U_y, U_z, mu, kappa, F)
         implicit none
         integer,       intent(in)  :: nEqn
         integer,       intent(in)  :: nGradEqn
         integer         , intent(in)  :: N(3)
         real(kind=RP),    intent(in)  :: Q  (1:nEqn, 0:N(1), 0:N(2), 0:N(3))
         real(kind=RP),    intent(in)  :: U_x(1:nGradEqn, 0:N(1), 0:N(2), 0:N(3) )
         real(kind=RP),    intent(in)  :: U_y(1:nGradEqn, 0:N(1), 0:N(2), 0:N(3) )
         real(kind=RP),    intent(in)  :: U_z(1:nGradEqn, 0:N(1), 0:N(2), 0:N(3) )
         real(kind=RP),    intent(in)  :: mu  (0:N(1), 0:N(2), 0:N(3))
         real(kind=RP),    intent(in)  :: kappa(0:N(1), 0:N(2), 0:N(3))
         real(kind=RP),    intent(out) :: F   (1:nEqn, 0:N(1), 0:N(2), 0:N(3), 1:NDIM )
!
!        ---------------
!        Local variables
!        ---------------
!
         integer       :: i , j , k

         do k = 0, N(3) ; do j = 0, N(2) ; do i = 0, N(1)
            F(INSRHO ,i,j,k,IX) = 0.0_RP
            F(INSU,i,j,k,IX)    = mu(i,j,k) * U_x(INSU,i,j,k)
            F(INSV,i,j,k,IX)    = mu(i,j,k) * U_x(INSV,i,j,k)
            F(INSW,i,j,k,IX)    = mu(i,j,k) * U_x(INSW,i,j,k) 
            F(INSP,i,j,k,IX)    = 0.0_RP
   
            F(INSRHO ,i,j,k,IY) = 0.0_RP
            F(INSU,i,j,k,IY)    = mu(i,j,k) * U_y(INSU,i,j,k)
            F(INSV,i,j,k,IY)    = mu(i,j,k) * U_y(INSV,i,j,k)
            F(INSW,i,j,k,IY)    = mu(i,j,k) * U_y(INSW,i,j,k) 
            F(INSP,i,j,k,IY)    = 0.0_RP
   
            F(INSRHO ,i,j,k,IZ) = 0.0_RP
            F(INSU,i,j,k,IZ)    = mu(i,j,k) * U_z(INSU,i,j,k)
            F(INSV,i,j,k,IZ)    = mu(i,j,k) * U_z(INSV,i,j,k)
            F(INSW,i,j,k,IZ)    = mu(i,j,k) * U_z(INSW,i,j,k) 
            F(INSP,i,j,k,IZ)    = 0.0_RP
   
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
      u = ABS( Q(INSU) )
      v = ABS( Q(INSV) )
      w = ABS( Q(INSW) )
      a = sqrt(max(max(u,v),w)**2 + 4.0_RP * thermodynamics % rho0c02/Q(INSRHO))
      
      eigen(1) = u + a
      eigen(2) = v + a
      eigen(3) = w + a
      
      END SUBROUTINE ComputeEigenvaluesForState
