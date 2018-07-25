!
!//////////////////////////////////////////////////////
!
!   @File:    TimeIntegratorDefinitions.f90
!   @Author:  Juan Manzanero (juan.manzanero@upm.es)
!   @Created: Tue Feb 13 14:26:34 2018
!   @Last revision date: Tue Jul  3 17:26:49 2018
!   @Last revision author: Juan Manzanero (juan.manzanero@upm.es)
!   @Last revision commit: 96905b05f7c99a4dc1a38da8202804d6dfef8cb3
!
!//////////////////////////////////////////////////////
!
module TimeIntegratorDefinitions
   use SMConstants
   implicit none

   private
   public   TimeStep_FCN

   abstract interface
      subroutine TimeStep_FCN( mesh, particles, t, BCFunctions, deltaT, ComputeTimeDerivative )
         use SMConstants
         use HexMeshClass
         use DGSEMClass
#if defined(NAVIERSTOKES) || defined(INCNS)
         use ParticlesClass, only: Particles_t
#endif
         IMPLICIT NONE
         type(HexMesh)              :: mesh
#if defined(NAVIERSTOKES) || defined(INCNS)
         type(Particles_t)          :: particles
#else
         logical                    :: particles
#endif
         REAL(KIND=RP)              :: t, deltaT
         type(BCFunctions_t), intent(in)  :: BCFunctions(no_of_BCsets)
         procedure(ComputeTimeDerivative_f) :: ComputeTimeDerivative
      end subroutine TimeStep_FCN
   end interface

end module TimeIntegratorDefinitions
