!
!//////////////////////////////////////////////////////
!
!   @File:    ArtificialViscosity.f90
!   @Author:  Juan Manzanero (juan.manzanero@upm.es)
!   @Created: Tue Aug 21 18:27:46 2018
!   @Last revision date:
!   @Last revision author:
!   @Last revision commit:
!
!//////////////////////////////////////////////////////
!
#include "Includes.h"
module ArtificialViscosity
   use SMConstants
   use PhysicsStorage_NS
   use VariableConversion_NS, only: Pressure
   use FluidData_NS, only: thermodynamics


   private
   public ComputeShockSensor

   real(kind=RP), parameter :: sBeta0 = 0.01_RP
   real(kind=RP), parameter :: eps_omega = 1.0e-10_RP

   contains
      subroutine ComputeShockSensor(Q, U_x, U_y, U_z, h, sBeta)
         implicit none
         real(kind=RP), intent(in)  :: Q(NCONS)
         real(kind=RP), intent(in)  :: U_x(NGRAD)
         real(kind=RP), intent(in)  :: U_y(NGRAD)
         real(kind=RP), intent(in)  :: U_z(NGRAD)
         real(kind=RP), intent(in)  :: h 
         real(kind=RP), intent(out) :: sBeta
!
!        ---------------
!        Local variables
!        ---------------
!
         real(kind=RP) :: sBetaMax, sTheta, sOmega, c, p, divV, rotV(3)

         sBetaMax = 2.0_RP / sqrt(POW2(thermodynamics % gamma) - 1.0_RP)

         p = Pressure(Q)
         c = sqrt(thermodynamics % gamma * p / Q(IRHO))

         divV = U_x(IGU) + U_y(IGV) + U_z(IGW)

         rotV(1) = U_y(IGW) - U_z(IGV)
         rotV(2) = U_z(IGU) - U_x(IGW)
         rotV(3) = U_x(IGV) - U_y(IGU)

         sTheta = -h * divV / c
         sOmega = (divV*divV)/(divV*divV + sum(rotV*rotV) + eps_omega)

         sBeta = sTheta * sOmega

         call SmoothedLimiter(sBeta, sBeta0, sBetaMax)

         sBeta = sBeta * sqrt(sum(POW2(Q(IRHOU:IRHOW))/(POW2(Q(IRHO)))) + c*c)
         
      end subroutine ComputeShockSensor

      subroutine SmoothedLimiter(s, s0, sMax)
         implicit none
         real(kind=RP), intent(inout) :: s
         real(kind=RP), intent(in)    :: s0
         real(kind=RP), intent(in)    :: sMax

         s = lmin(lmax(s-s0) - sMax) + sMax

      end subroutine SmoothedLimiter

      pure function Lmax(s)
         implicit none
         real(kind=RP), intent(in) :: s
         real(kind=RP)             :: Lmax
         real(kind=RP), parameter  :: b = 100.0_RP

         Lmax = (s*atan(b*s) - atan(b))/PI + 0.5_RP*(s+1)

      end function

      pure function Lmin(s)
         implicit none
         real(kind=RP), intent(in) :: s
         real(kind=RP)             :: Lmin

         Lmin = s - Lmax(s)
      
      end function


end module ArtificialViscosity
