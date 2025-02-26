#include "Includes.h"
module VariableConversion_SLR_INS_V04
    use SMConstants
    use PhysicsStorage_SLR_INS_V04
    use FluidData_SLR_INS_V04
    implicit none
 
    private
    public   SLR_INS_V04GradientVariables
    public   GetViscosity

    contains
!
! /////////////////////////////////////////////////////////////////////
!
!---------------------------------------------------------------------
!! GradientValuesForQ takes the solution (Q) values and returns the
!! quantities of which the gradients will be taken.
!---------------------------------------------------------------------
!
    pure subroutine SLR_INS_V04GradientVariables(nEqn, nGrad, Q, U, rho_)
        implicit none
        integer, intent(in)        :: nEqn, nGrad
        real(kind=RP), intent(in)  :: Q(nEqn)
        real(kind=RP), intent(out) :: U(nGrad)
        real(kind=RP), intent(in), optional :: rho_
        
        U = Q

    end subroutine SLR_INS_V04GradientVariables
    
    pure subroutine GetViscosity(phi, mu)
        implicit none
        real(kind=RP), intent(in)   :: phi
        real(kind=RP), intent(out)  :: mu

        !//#
        !mu = dimensionless % mu
        mu =1.0

    end subroutine GetViscosity

end module VariableConversion_SLR_INS_V04