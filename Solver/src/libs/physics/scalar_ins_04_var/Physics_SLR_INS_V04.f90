#include "Includes.h"
!  **************
   module Physics_SLR_INS_V04
!  **************
!
      use SMConstants
      use PhysicsStorage_SLR_INS_V04
      use VariableConversion_SLR_INS_V04
      use FluidData_SLR_INS_V04
      implicit none
      
      private
      public  slr_ins_v04ViscousFlux
      public iEulerFlux_INS
      public iEulerFlux_INSTest
      public iEulerFlux_INS_burges

      

      ! procedure :: iEulerFlux_INS                      => iEulerFlux_INSTest

!
!     ========
      CONTAINS 
!     ========
! !
      pure subroutine iEulerFlux_INS(Q, F)
      implicit none
      real(kind=RP), intent(in)   :: Q(1:N_INS)
      real(kind=RP), intent(out)  :: F(1:N_INS, 1:NDIM)
      ! real(kind=RP), intent(in), optional :: rho_
!
!        ---------------
!        Local variables
!        ---------------
!
            real(kind=RP) :: rho, invRho
!     
!              X-Flux
!              ------         
            F(1, IX) = Q(1)*Q(1)
            F(2, IX) = Q(1)*Q(2)
            F(3, IX) = Q(1)*Q(3)

!              Y-Flux
! !              ------
            F(1, IY) = Q(2)*Q(1)
            F(2, IY) = Q(2)*Q(2)
            F(3, IY) = Q(2)*Q(3)
!     
!              Z-Flux
!              ------
            F(1,IZ) = Q(3)*Q(1)
            F(2,IZ) = Q(3)*Q(2)
            F(3,IZ) = Q(3)*Q(3)
   
   end subroutine iEulerFlux_INS


   pure subroutine iEulerFlux_INS_burges(Q, F)
      implicit none
      real(kind=RP), intent(in)   :: Q(1:N_INS)
      real(kind=RP), intent(out)  :: F(1:N_INS, 1:NDIM)
      ! real(kind=RP), intent(in), optional :: rho_
!
!        ---------------
!        Local variables
!        ---------------
!
      real(kind=RP) :: rho, invRho
!              X-Flux
!              ------         
      F(1, IX) = Q(1)*Q(1)
      F(2, IX) = Q(1)*Q(2)
      F(3, IX) = Q(1)*Q(3)
!              Y-Flux
! !              ------
      F(1, IY) = Q(2)*Q(1)
      F(2, IY) = Q(2)*Q(2)
      F(3, IY) = Q(2)*Q(3)
!              Z-Flux
!              ------
      F(1,IZ) = Q(3)*Q(1)
      F(2,IZ) = Q(3)*Q(2)
      F(3,IZ) = Q(3)*Q(3)

   end subroutine iEulerFlux_INS_burges

   pure subroutine iEulerFlux_INSTest(Q, F)
      implicit none
      real(kind=RP), intent(in)   :: Q(1:N_INS)
      real(kind=RP), intent(out)  :: F(1:N_INS, 1:NDIM)
      ! real(kind=RP), intent(in), optional :: rho_
!
!        ---------------
!        Local variables
!        ---------------
!
            real(kind=RP) :: rho, invRho

!              X-Flux
!              ------         
            F(1, IX) = Q(1)*0.5_RP
            F(2, IX) = Q(2)*0.5_RP
            F(3, IX) = Q(3)*0.5_RP

            !              Y-Flux
! !              ------
            F(1, IY) = Q(1)*0.0_RP
            F(2, IY) = Q(2)*0.0_RP
            F(3, IY) = Q(3)*0.0_RP
!     
!              Z-Flux
!              ------
            F(1,IZ) = Q(1)*0.5_RP
            F(2,IZ) = Q(2)*0.5_RP
            F(3,IZ) = Q(3)*0.5_RP
   
   end subroutine iEulerFlux_INSTest
!
!//////////////////////////////////////////////////////////////////////////////////////////
!
!        VISCOUS FLUXES
!         --------------
!
!//////////////////////////////////////////////////////////////////////////////////////////
!
      pure subroutine slr_ins_v04ViscousFlux(nEqn, nGradEqn, Q, U_x, U_y, U_z, mu, beta, kappa, F)
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

         ! Change the indecies to slr_ins_v04
         F(1,IX) = mu * U_x(1)
         F(2,IX) = mu * U_x(2)
         F(3,IX) = mu * U_x(3)
         F(4,IX) = mu * U_x(4)

         F(1,IY)  = mu * U_y(1)
         F(2,IY)  = mu * U_y(2)
         F(3,IY)  = mu * U_y(3)
         F(4,IY)  = mu * U_y(4)

         F(1,IZ)  = mu * U_z(1)
         F(2,IZ)  = mu * U_z(2)
         F(3,IZ)  = mu * U_z(3)
         F(4,IZ)  = mu * U_z(4)


      end subroutine slr_ins_v04ViscousFlux
 
   END Module Physics_SLR_INS_V04
