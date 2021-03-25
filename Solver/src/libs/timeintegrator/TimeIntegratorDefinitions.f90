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
   public   TimeStep_FCN, ComputePseudoTimeDerivative

   abstract interface
      subroutine TimeStep_FCN( mesh, particles, t, deltaT, ComputeTimeDerivative , dt_vec, dts, global_dt )
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
         ! Optional arguments:
         real(kind=RP), allocatable, dimension(:), intent(in), optional :: dt_vec ! dt array for Local Time Stepping preconditioning
         ! Dual (pseudo) time stepping arguments (also optional):
         logical, intent(in), optional :: dts 
         real(kind=RP), intent(in), optional :: global_dt 
      end subroutine TimeStep_FCN
   end interface
!========
   contains
!========
!
!////////////////////////////////////////////////////////////////////////
!
   subroutine ComputePseudoTimeDerivative(mesh, tk, global_dt)
      use SMConstants
      use HexMeshClass
!-----Arguments-----------------------------------------------------------
      type(HexMesh)                                        :: mesh
      real(KIND=RP), intent(in)                            :: tk
      real(kind=RP), intent(in)                            :: global_dt 
!-----Local-Variables-----------------------------------------------------
      integer       :: id
!-------------------------------------------------------------------------

!$omp parallel do schedule(runtime)
   do id = 1, SIZE( mesh % elements )
      mesh % elements(id) % storage % Qdot = mesh % elements(id) % storage % Qdot - &
         (mesh % elements(id) % storage % Q - mesh % elements(id) % storage % prevQ(1) % Q) / global_dt
   end do ! id
!$omp end parallel do

   end subroutine
!
!////////////////////////////////////////////////////////////////////////
!
end module TimeIntegratorDefinitions
