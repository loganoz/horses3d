!
!//////////////////////////////////////////////////////
!
!   @File:    VariableConversion.f90
!   @Author:  Juan Manzanero (juan.manzanero@upm.es)
!   @Created: Tue Jan 16 11:59:44 2018
!   @Last revision date:
!   @Last revision author:
!   @Last revision commit:
!
!//////////////////////////////////////////////////////
!
#include "Includes.h"
module VariableConversion
   use SMConstants
   use PhysicsStorage
   implicit none

   private
   public   GradientValuesForQ

   interface GradientValuesForQ
       module procedure GradientValuesForQ_0D , GradientValuesForQ_3D
   end interface GradientValuesForQ

   contains
      pure SUBROUTINE GradientValuesForQ_0D( Q, U )
      IMPLICIT NONE
!
!     ---------
!     Arguments
!     ---------
!
      REAL(KIND=RP), DIMENSION(N_EQN)     , INTENT(IN)  :: Q
      REAL(KIND=RP), DIMENSION(N_GRAD_EQN), INTENT(OUT) :: U

      U = Q

      END SUBROUTINE GradientValuesForQ_0D

      pure SUBROUTINE GradientValuesForQ_3D( Nx, Ny, Nz, Q, U )
      IMPLICIT NONE
!
!     ---------
!     Arguments
!     ---------
!
      integer,       intent(in)  :: Nx, Ny, Nz
      REAL(KIND=RP), INTENT(IN)  :: Q(1:NCONS, 0:Nx, 0:Ny, 0:Nz)
      REAL(KIND=RP), INTENT(OUT) :: U(1:N_GRAD_EQN, 0:Nx, 0:Ny, 0:Nz)

      U = Q

      END SUBROUTINE GradientValuesForQ_3D

end module VariableConversion
