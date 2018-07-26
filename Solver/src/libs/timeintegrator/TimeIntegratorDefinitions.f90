!
!//////////////////////////////////////////////////////
!
!   @File:    TimeIntegratorDefinitions.f90
!   @Author:  Juan Manzanero (juan.manzanero@upm.es)
!   @Created: Tue Feb 13 14:26:34 2018
!   @Last revision date: Wed Jul 25 17:15:49 2018
!   @Last revision author: Juan Manzanero (juan.manzanero@upm.es)
!   @Last revision commit: d886ff7a7d37081df645692157131f3ecc98f761
!
!//////////////////////////////////////////////////////
!
module TimeIntegratorDefinitions
   use SMConstants
   implicit none

   private
   public   TimeStep_FCN

   abstract interface
      subroutine TimeStep_FCN( mesh, particles, t, deltaT, ComputeTimeDerivative )
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
         procedure(ComputeTimeDerivative_f) :: ComputeTimeDerivative
      end subroutine TimeStep_FCN
   end interface

end module TimeIntegratorDefinitions
