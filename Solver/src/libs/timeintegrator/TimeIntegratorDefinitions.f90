!
!//////////////////////////////////////////////////////
!
!   @File:    TimeIntegratorDefinitions.f90
!   @Author:  Juan Manzanero (juan.manzanero@upm.es)
!   @Created: Tue Feb 13 14:26:34 2018
!   @Last revision date: Tue Apr 10 12:12:12 2018
!   @Last revision author: Juan (juan.manzanero@upm.es)
!   @Last revision commit: 00cf81fc8040d50de254463075cfb5043c152064
!
!//////////////////////////////////////////////////////
!
module TimeIntegratorDefinitions
   use SMConstants
   implicit none

   private
   public   TimeStep_FCN

   abstract interface
      subroutine TimeStep_FCN( mesh, t, BCFunctions, deltaT, ComputeTimeDerivative )
         use SMConstants
         use HexMeshClass
         use DGSEMClass
         IMPLICIT NONE
         type(HexMesh)              :: mesh
         REAL(KIND=RP)              :: t, deltaT
         type(BCFunctions_t), intent(in)  :: BCFunctions(no_of_BCsets)
         procedure(ComputeQDot_FCN) :: ComputeTimeDerivative
      end subroutine TimeStep_FCN
   end interface

end module TimeIntegratorDefinitions
