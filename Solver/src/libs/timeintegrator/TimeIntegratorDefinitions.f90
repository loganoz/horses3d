!
!//////////////////////////////////////////////////////
!
!   @File:    TimeIntegratorDefinitions.f90
!   @Author:  Juan Manzanero (juan.manzanero@upm.es)
!   @Created: Tue Feb 13 14:26:34 2018
!   @Last revision date: Tue Feb 13 19:37:44 2018
!   @Last revision author: Juan (juan.manzanero@upm.es)
!   @Last revision commit: c01958bbb74b2de9252027cd1c501fe081a58ef2
!
!//////////////////////////////////////////////////////
!
module TimeIntegratorDefinitions
   use SMConstants
   implicit none

   private
   public   TimeStep_FCN

   abstract interface
      subroutine TimeStep_FCN( mesh, t, externalState, externalGradients, deltaT, ComputeTimeDerivative )
         use SMConstants
         use HexMeshClass
         use DGSEMClass
         IMPLICIT NONE
         type(HexMesh)              :: mesh
         REAL(KIND=RP)              :: t, deltaT
         external                   :: externalState, externalGradients
         procedure(ComputeQDot_FCN) :: ComputeTimeDerivative
      end subroutine TimeStep_FCN
   end interface

end module TimeIntegratorDefinitions
