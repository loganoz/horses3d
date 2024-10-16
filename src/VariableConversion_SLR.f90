#include "Includes.h"
module VariableConversion_SLR
    use SMConstants
    use PhysicsStorage_SLR
    use FluidData_SLR
    implicit none
 
    private
    public   SLRGradientVariables
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
    pure subroutine SLRGradientVariables(nEqn, nGrad, Q, U, rho_)
        implicit none
        integer, intent(in)        :: nEqn, nGrad
        real(kind=RP), intent(in)  :: Q(nEqn)
        real(kind=RP), intent(out) :: U(nGrad)
        real(kind=RP), intent(in), optional :: rho_
        
        U = Q

    end subroutine SLRGradientVariables
    
    pure subroutine GetViscosity(phi, mu)
        implicit none
        real(kind=RP), intent(in)   :: phi
        real(kind=RP), intent(out)  :: mu

        !//#
        !mu = dimensionless % mu
        mu =1.0

    end subroutine GetViscosity

end module VariableConversion_SLR