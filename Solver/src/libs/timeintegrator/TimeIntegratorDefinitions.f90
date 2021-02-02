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
#include "Includes.h"
module TimeIntegratorDefinitions
   use SMConstants
   implicit none

   private
   public   TimeStep_FCN

   abstract interface
      subroutine TimeStep_FCN( mesh, particles, t, deltaT, ComputeTimeDerivative , dt_vec)
         use SMConstants
         use HexMeshClass
         use DGSEMClass
#ifdef FLOW
         use ParticlesClass, only: Particles_t
#endif
         IMPLICIT NONE
         type(HexMesh)              :: mesh
#ifdef FLOW
         type(Particles_t)          :: particles
#else
         logical                    :: particles
#endif
         REAL(KIND=RP)              :: t, deltaT
         procedure(ComputeTimeDerivative_f) :: ComputeTimeDerivative
         real(kind=RP), allocatable, dimension(:), intent(in), optional :: dt_vec
      end subroutine TimeStep_FCN
   end interface

end module TimeIntegratorDefinitions
