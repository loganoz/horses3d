#include "Includes.h"
module TimeIntegratorDefinitions
   use SMConstants
   implicit none

   private
   public   TimeStep_FCN, ComputePseudoTimeDerivative, bdf_order, NOTDEF_JACOBIAN, NUMERICAL_JACOBIAN, ANALYTICAL_JACOBIAN
   public JACOBIAN_MATRIX_NONE, JACOBIAN_MATRIX_DENSE, JACOBIAN_MATRIX_CSR

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

   integer       :: bdf_order       ! BDF order specified by user
   integer, parameter :: NOTDEF_JACOBIAN     = 0
   integer, parameter :: NUMERICAL_JACOBIAN  = 1
   integer, parameter :: ANALYTICAL_JACOBIAN = 2

   integer, parameter :: JACOBIAN_MATRIX_NONE  = 0
   integer, parameter :: JACOBIAN_MATRIX_DENSE = 1
   integer, parameter :: JACOBIAN_MATRIX_CSR   = 2
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
      ! integer,       intent(in)                            :: bdf_order 
!-----Local-Variables-----------------------------------------------------
      integer       :: id, i
      real(kind=RP) :: invdt 
      real(kind=RP), parameter :: BDFCoeff(6,5) = &
      !                    a_1             a_2     a_3           a_4             a_5          a_6
               reshape( (/ 1.0_RP        , -1._RP, 0._RP       , 0._RP         , 0._RP      , 0._RP        ,  &   ! BDF1
                           1.5_RP        , -2._RP, 0.5_RP      , 0._RP         , 0._RP      , 0._RP        ,  &   ! BDF2
                           11._RP/6_RP   , -3._RP, 3._RP/2._RP , -1._RP/3._RP  , 0._RP      , 0._RP        ,  &   ! BDF3
                           25._RP/12_RP  , -4._RP, 3._RP       , -4._RP/3._RP  , 1._RP/4._RP, 0._RP        ,  &   ! BDF4
                           137._RP/60_RP , -5._RP, 5._RP       , -10._RP/3._RP , 5._RP/4._RP, -1._RP/5._RP /) &   ! BDF5
                                                                                                            , (/6,5/) )
!-------------------------------------------------------------------------

      invdt = 1._RP/global_dt

      select case (bdf_order)
      case (1)
!$omp parallel do schedule(runtime)
         do id = 1, SIZE( mesh % elements )
            mesh % elements(id) % storage % Qdot = mesh % elements(id) % storage % Qdot - &
               (BDFCoeff(1, bdf_order) * mesh % elements(id) % storage % Q + &
                BDFCoeff(2, bdf_order) * mesh % elements(id) % storage % prevQ(1) % Q) * invdt
         end do ! id
!$omp end parallel do
      case (2)
!$omp parallel do schedule(runtime)
         do id = 1, SIZE( mesh % elements )
            mesh % elements(id) % storage % Qdot = mesh % elements(id) % storage % Qdot - &
               (BDFCoeff(1, bdf_order) * mesh % elements(id) % storage % Q + &
                BDFCoeff(2, bdf_order) * mesh % elements(id) % storage % prevQ(2) % Q + &
                BDFCoeff(3, bdf_order) * mesh % elements(id) % storage % prevQ(1) % Q) * invdt
         end do ! id
!$omp end parallel do
      case (3)
!$omp parallel do schedule(runtime)
         do id = 1, SIZE( mesh % elements )
            mesh % elements(id) % storage % Qdot = mesh % elements(id) % storage % Qdot - &
               (BDFCoeff(1, bdf_order) * mesh % elements(id) % storage % Q + &
                BDFCoeff(2, bdf_order) * mesh % elements(id) % storage % prevQ(3) % Q + &
                BDFCoeff(3, bdf_order) * mesh % elements(id) % storage % prevQ(2) % Q + &
                BDFCoeff(4, bdf_order) * mesh % elements(id) % storage % prevQ(1) % Q) * invdt
         end do ! id
!$omp end parallel do
      case (4)
!$omp parallel do schedule(runtime)
         do id = 1, SIZE( mesh % elements )
            mesh % elements(id) % storage % Qdot = mesh % elements(id) % storage % Qdot - &
               (BDFCoeff(1, bdf_order) * mesh % elements(id) % storage % Q + &
                BDFCoeff(2, bdf_order) * mesh % elements(id) % storage % prevQ(4) % Q + &
                BDFCoeff(3, bdf_order) * mesh % elements(id) % storage % prevQ(3) % Q + &
                BDFCoeff(4, bdf_order) * mesh % elements(id) % storage % prevQ(2) % Q + &
                BDFCoeff(5, bdf_order) * mesh % elements(id) % storage % prevQ(1) % Q) * invdt
         end do ! id
!$omp end parallel do
      case (5)
!$omp parallel do schedule(runtime)
         do id = 1, SIZE( mesh % elements )
            mesh % elements(id) % storage % Qdot = mesh % elements(id) % storage % Qdot - &
               (BDFCoeff(1, bdf_order) * mesh % elements(id) % storage % Q + &
                BDFCoeff(2, bdf_order) * mesh % elements(id) % storage % prevQ(5) % Q + &
                BDFCoeff(3, bdf_order) * mesh % elements(id) % storage % prevQ(4) % Q + &
                BDFCoeff(4, bdf_order) * mesh % elements(id) % storage % prevQ(3) % Q + &
                BDFCoeff(5, bdf_order) * mesh % elements(id) % storage % prevQ(2) % Q + &
                BDFCoeff(6, bdf_order) * mesh % elements(id) % storage % prevQ(1) % Q) * invdt
         end do ! id
!$omp end parallel do
      end select


   end subroutine
!
!////////////////////////////////////////////////////////////////////////
!
end module TimeIntegratorDefinitions