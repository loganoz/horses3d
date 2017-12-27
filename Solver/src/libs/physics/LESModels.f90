!
!//////////////////////////////////////////////////////
!
!   @File:    LESModels.f90
!   @Author:  Juan Manzanero (juan.manzanero@upm.es)
!   @Created: Wed Dec 27 17:44:13 2017
!   @Last revision date:
!   @Last revision author:
!   @Last revision commit:
!
!//////////////////////////////////////////////////////
!
#include "Includes.h"
module LESModels
   use SMConstants
   use PhysicsStorage
   implicit none

   private
   public BasicSmagorinskySGSTensor

   real(kind=RP), parameter   :: CS = 0.2_RP

   interface BasicSmagorinskySGSTensor
      module procedure BasicSmagorinskySGSTensor3D
      module procedure BasicSmagorinskySGSTensor0D
   end interface BasicSmagorinskySGSTensor   

   
   contains
      subroutine BasicSmagorinskySGSTensor3D(delta, N, U_x, U_y, U_z, tau)
         implicit none
         real(kind=RP), intent(in)     :: delta
         integer,       intent(in)     :: N(3)
         real(kind=RP), intent(in)     :: U_x(NGRAD, 0:N(1), 0:N(2), 0:N(3))
         real(kind=RP), intent(in)     :: U_y(NGRAD, 0:N(1), 0:N(2), 0:N(3))
         real(kind=RP), intent(in)     :: U_z(NGRAD, 0:N(1), 0:N(2), 0:N(3))
         real(kind=RP), intent(out)    :: tau(NDIM, NDIM, 0:N(1), 0:N(2), 0:N(3))
!
!        ---------------
!        Local variables
!        ---------------
!
         integer     :: i, j, k
         real(kind=RP)  :: S(NDIM, NDIM)
         real(kind=RP)  :: normS, divV

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
            normS = sqrt( 2.0_RP * product(S*S) )
!
!           Remove the volumetric deformation tensor
!           ----------------------------------------
            S(1,1) = S(1,1) - 1.0_RP / 3.0_RP * divV
            S(2,2) = S(2,2) - 1.0_RP / 3.0_RP * divV
            S(3,3) = S(3,3) - 1.0_RP / 3.0_RP * divV
!
!           Compute the SGS tensor
!           ----------------------
            tau(:,:,i,j,k) = -2.0_RP * POW2(CS*delta) * normS * S

         end do         ; end do          ; end do
         
      end subroutine BasicSmagorinskySGSTensor3D

      subroutine BasicSmagorinskySGSTensor0D(delta, U_x, U_y, U_z, tau)
         implicit none
         real(kind=RP), intent(in)     :: delta
         real(kind=RP), intent(in)     :: U_x(NGRAD)
         real(kind=RP), intent(in)     :: U_y(NGRAD)
         real(kind=RP), intent(in)     :: U_z(NGRAD)
         real(kind=RP), intent(out)    :: tau(NDIM, NDIM)
!
!        ---------------
!        Local variables
!        ---------------
!
         real(kind=RP)  :: S(NDIM, NDIM)
         real(kind=RP)  :: normS, divV
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
         normS = sqrt( 2.0_RP * product(S*S) )
!
!        Remove the volumetric deformation tensor
!        ----------------------------------------
         S(1,1) = S(1,1) - 1.0_RP / 3.0_RP * divV
         S(2,2) = S(2,2) - 1.0_RP / 3.0_RP * divV
         S(3,3) = S(3,3) - 1.0_RP / 3.0_RP * divV
!
!        Compute the SGS tensor
!        ----------------------
         tau = -2.0_RP * POW2(CS*delta) * normS * S

      end subroutine BasicSmagorinskySGSTensor0D

end module LESModels
