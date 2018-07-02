!
!//////////////////////////////////////////////////////
!
!   @File:    TimeIntegratorDefinitions.f90
!   @Author:  Juan Manzanero (juan.manzanero@upm.es)
!   @Created: Tue Feb 13 14:26:34 2018
!   @Last revision date: Sat Jun 23 10:20:39 2018
!   @Last revision author: Juan Manzanero (juan.manzanero@upm.es)
!   @Last revision commit: fce351220409e80ce5df1949249c2b870dd847aa
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
         procedure(ComputeQDot_FCN) :: ComputeTimeDerivative
      end subroutine TimeStep_FCN
   end interface

end module TimeIntegratorDefinitions
