!
!//////////////////////////////////////////////////////
!
!   @File:    RiemannSolvers_NS.f90
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
   module RiemannSolvers_NS
      use SMConstants
      use Physics_NS
      use PhysicsStorage_NS
      use FluidData, only: equationOfState, getThermalConductivity
      contains
         subroutine SetRiemannSolver(which, splitType)
            integer, intent(in)  :: which
            integer, intent(in)  :: splitType
         end subroutine SetRiemannSolver
!
!     ////////////////////////////////////////////////////////////////////////////////////////
!
      SUBROUTINE RiemannSolver( QLeft, QRight, nHat, t1, t2, flux )
         IMPLICIT NONE
!
!        ---------
!        Arguments
!        ---------
!
         REAL(KIND=RP), DIMENSION(NCONS)  :: Qleft, Qright, flux
         REAL(KIND=RP), DIMENSION(3)      :: nHat, t1, t2
         
         flux = 0.5_RP*(Qleft + Qright)*( nHat(1) + nHat(2) + nHat(3) )
      
      END SUBROUTINE RiemannSolver
      
      SUBROUTINE RiemannSolver_dFdQ(ql,qr,nHat,dfdq_num,side)
         IMPLICIT NONE
!
!        ---------
!        Arguments
!        ---------
!
         REAL(KIND=RP), DIMENSION(NCONS)        :: Ql, Qr
         REAL(KIND=RP), DIMENSION(NCONS,NCONS)  :: dfdq_num
         integer                                :: side
         REAL(KIND=RP), DIMENSION(3)            :: nHat
         
      
      END SUBROUTINE RiemannSolver_dFdQ
   end module RiemannSolvers_NS

