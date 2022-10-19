!
!//////////////////////////////////////////////////////
!
!   @File:    VariableConversion_NS.f90
!   @Author:  Juan Manzanero (juan.manzanero@upm.es)
!   @Created: Mon May 14 19:03:30 2018
!   @Last revision date:
!   @Last revision author:
!   @Last revision commit:
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
module VariableConversion_NS
   use SMConstants
   use PhysicsStorage_NS
   use FluidData_NS
   implicit none

   private
   public   Pressure, Temperature, NSGradientValuesForQ
   public   NSGradientValuesForQ_0D, NSGradientValuesForQ_3D
   public   getPrimitiveVariables, getEntropyVariables
   public   getRoeVariables

   interface NSGradientValuesForQ
       module procedure NSGradientValuesForQ_0D , NSGradientValuesForQ_3D
   end interface NSGradientValuesForQ

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
      
      P = thermodynamics % gammaMinus1*(Q(5) - 0.5_RP*(Q(2)**2 + Q(3)**2 + Q(4)**2)/Q(1))

      end function Pressure
!
!
! /////////////////////////////////////////////////////////////////////
!
!---------------------------------------------------------------------
!! Compute the temperature from the state variables
!---------------------------------------------------------------------
!
      PURE function Temperature(Q) RESULT(T)
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
      REAL(KIND=RP) :: T
!
      T = dimensionless % gammaM2*Pressure(Q)/Q(1)

      end function Temperature
!
! /////////////////////////////////////////////////////////////////////
!
!---------------------------------------------------------------------
!! GradientValuesForQ takes the solution (Q) values and returns the
!! quantities of which the gradients will be taken.
!---------------------------------------------------------------------
!
      pure subroutine NSGradientValuesForQ_0D( nEqn, nGrad, Q, U )
         implicit none
         integer, intent(in)        :: nEqn, nGrad
         real(kind=RP), intent(in)  :: Q(nEqn)
         real(kind=RP), intent(out) :: U(nGrad)
!
!        ---------------
!        Local Variables
!        ---------------
!     

      U(1) = Q(1)
      U(2) = Q(2)
      U(3) = Q(3)
      U(4) = Q(4)

      end subroutine NSGradientValuesForQ_0D

      pure subroutine NSGradientValuesForQ_3D( nEqn, nGrad, Nx, Ny, Nz, Q, U )
         implicit none
         integer,       intent(in)  :: nEqn, nGrad, Nx, Ny, Nz
         real(kind=RP), intent(in)  :: Q(1:nEqn,  0:Nx, 0:Ny, 0:Nz)
         real(kind=RP), intent(out) :: U(1:nGrad, 0:Nx, 0:Ny, 0:Nz)
!
!        ---------------
!        Local Variables
!        ---------------
!     
         integer     :: i, j, k

         associate ( gammaM2 => dimensionless % gammaM2, &
                     gammaMinus1 => thermodynamics % gammaMinus1 ) 
      

      U(1,:,:,:) = Q(1,:,:,:)
      U(2,:,:,:) = Q(2,:,:,:)
      U(3,:,:,:) = Q(3,:,:,:)
      U(4,:,:,:) = Q(4,:,:,:)
         end associate

      end subroutine NSGradientValuesForQ_3D
!
! /////////////////////////////////////////////////////////////////////
!
      pure subroutine getPrimitiveVariables(U,V)
!
!        **************************************
!           Primitive variables are:
!              V = [invRho, u,v,w,p,T,a^2]
!        **************************************
!
         implicit none
         real(kind=RP), intent(in)  :: U(NCONS)
         real(kind=RP), intent(out) :: V(NPRIM)
!
!        ---------------
!        Local variables
!        ---------------
!
         real(kind=RP)  :: invRho

         invRho = 1.0_RP / U(IRHO)

         V(IPIRHO) = invRho
         V(IPU) = U(IRHOU) * invRho
         V(IPV) = U(IRHOV) * invRho
         V(IPW) = U(IRHOW) * invRho
         V(IPP) = thermodynamics % gammaMinus1 * ( U(IRHOE) &
                  - 0.5_RP * (V(IPU)*U(IRHOU) + V(IPV)*U(IRHOV) + V(IPW)*U(IRHOW)))
         V(IPT) = V(IPP) * dimensionless % gammaM2 * invRho 
         V(IPA2) = thermodynamics % gamma * V(IPP) * invRho

      end subroutine getPrimitiveVariables

      pure subroutine getEntropyVariables(U,p,invRho,S)
         implicit none
         real(kind=RP), intent(in)  :: U(NCONS)
         real(kind=RP), intent(in)  :: p, invRho
         real(kind=RP), intent(out) :: S(NCONS)
!
!        ---------------
!        Local variables
!        ---------------
!
         real(kind=RP)  :: invP
         real(kind=RP)  :: entropy

         invP = 1.0_RP / p

         entropy = log(p * (invRho ** thermodynamics % gamma))
        
         S(1) =   (thermodynamics % gamma - entropy) / (thermodynamics % gammaMinus1) &
                - 0.5_RP * invRho * (POW2(U(IRHOU))+POW2(U(IRHOV))+POW2(U(IRHOW))) * invP

         S(2) = U(IRHOU) * invP
         S(3) = U(IRHOV) * invP
         S(4) = U(IRHOW) * invP
         S(5) = -U(IRHO) * invP

      end subroutine getEntropyVariables

      pure subroutine getRoeVariables(QL, QR, VL, VR, rho, u, v, w, V2, H, a)
!
!        ***************************************************
!           Roe variables are: [rho, u, v, w, H, a]
!        ***************************************************
!           
         implicit none
         real(kind=RP), intent(in)  :: QL(NCONS)
         real(kind=RP), intent(in)  :: QR(NCONS)
         real(kind=RP), intent(in)  :: VL(NPRIM)
         real(kind=RP), intent(in)  :: VR(NPRIM)
         real(kind=RP), intent(out) :: rho, u, v, w, V2, H, a
!
!        ---------------
!        Local variables
!        ---------------
!
         real(kind=RP)  :: HL, HR
         real(kind=RP)  :: sqrtRhoL, sqrtRhoR
         real(kind=RP)  :: invSumSqrtRhoLR

         associate(gamma => thermodynamics % gamma, &
                   gm1   => thermodynamics % gammaMinus1)

         sqrtRhoL = sqrt(QL(IRHO))  ; sqrtRhoR = sqrt(QR(IRHO))
         invSumSqrtRhoLR = 1.0_RP / (sqrtRhoL + sqrtRhoR)
!
!        Here the enthalpy is defined as rhoH = gogm1 p + 0.5 rho V^2 = p + rhoe
!        -----------------------------------------------------------------------
         HL = (VL(IPP) + QL(IRHOE)) * VL(IPIRHO)
         HR = (VR(IPP) + QR(IRHOE)) * VR(IPIRHO)

         rho = sqrtRhoL * sqrtRhoR
         u   = (sqrtRhoL * VL(IPU) + sqrtRhoR * VR(IPU))*invSumSqrtRhoLR
         v   = (sqrtRhoL * VL(IPV) + sqrtRhoR * VR(IPV))*invSumSqrtRhoLR
         w   = (sqrtRhoL * VL(IPW) + sqrtRhoR * VR(IPW))*invSumSqrtRhoLR
         H   = (sqrtRhoL * HL      + sqrtRhoR * HR     )*invSumSqrtRhoLR
         V2  = POW2(u) + POW2(v) + POW2(w)
         a   = sqrt(gm1*(H - 0.5_RP * V2))

         end associate

      end subroutine getRoeVariables

end module VariableConversion_NS
