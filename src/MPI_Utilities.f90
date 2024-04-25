#include "Includes.h"
module MPI_Utilities
   use SMConstants
   use MPI_Process_Info
#if defined(_HAS_MPI_)
   use mpi
#endif

   private
   public infNorm, L2Norm, MPI_SumAll, MPI_MinMax

   interface MPI_SumAll
      module procedure MPI_SumAll_int_r0
      module procedure MPI_SumAll_int_r1
      module procedure MPI_SumAll_int_r2
      module procedure MPI_SumAll_int_r3
      module procedure MPI_SumAll_real_r0
      module procedure MPI_SumAll_real_r1
      module procedure MPI_SumAll_real_r2
      module procedure MPI_SumAll_real_r3
   end interface

   interface MPI_MinMax
      module procedure MPI_MinMax_real_r0
      module procedure MPI_MinMax_real_r1
      module procedure MPI_MinMax_real_r2
      module procedure MPI_MinMax_real_r3
   end interface
!
!  ========
   contains
!  ========
!
   function L2Norm(x) result(norm)
      implicit none
      !-arguments---------------------------------------
      real(kind=RP), intent(in) :: x(:)
      real(kind=RP)             :: norm
      !-local-variables---------------------------------
      integer       :: ierr
      !-------------------------------------------------

      norm = sum(x**2)

#if defined(_HAS_MPI_)
      if (MPI_Process % doMPIAction) then
         call MPI_Allreduce(MPI_IN_PLACE, norm, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD, ierr)
      end if
#endif

      norm = sqrt(norm)

   end function L2Norm

   function infNorm(x) result(norm)
      implicit none
      !-arguments---------------------------------------
      real(kind=RP), intent(in) :: x(:)
      real(kind=RP)             :: norm
      !-local-variables---------------------------------
      integer       :: ierr
      !-------------------------------------------------

      norm = maxval(abs(x))

#if defined(_HAS_MPI_)
      if (MPI_Process % doMPIAction) then
         call MPI_Allreduce(MPI_IN_PLACE, norm, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD, ierr)
      end if
#endif

   end function infNorm

   subroutine MPI_SumAll_int_r0(x)
      implicit none
      !-arguments---------------------------------------
      integer, intent(inout) :: x
      !-local-variables---------------------------------
      integer :: ierr
      !-------------------------------------------------

#if defined(_HAS_MPI_)
      call MPI_Allreduce(MPI_IN_PLACE, x, 1, MPI_INTEGER, MPI_SUM, MPI_COMM_WORLD, ierr)
#endif
   end subroutine MPI_SumAll_int_r0

   subroutine MPI_SumAll_int_r1(x)
      implicit none
      !-arguments---------------------------------------
      integer, intent(inout) :: x(:)
      !-local-variables---------------------------------
      integer :: ierr
      !-------------------------------------------------

#if defined(_HAS_MPI_)
      call MPI_Allreduce(MPI_IN_PLACE, x, size(x), MPI_INTEGER, MPI_SUM, MPI_COMM_WORLD, ierr)
#endif
   end subroutine MPI_SumAll_int_r1

   subroutine MPI_SumAll_int_r2(x)
      implicit none
      !-arguments---------------------------------------
      integer, intent(inout) :: x(:,:)
      !-local-variables---------------------------------
      integer :: ierr
      !-------------------------------------------------

#if defined(_HAS_MPI_)
      call MPI_Allreduce(MPI_IN_PLACE, x, size(x), MPI_INTEGER, MPI_SUM, MPI_COMM_WORLD, ierr)
#endif
   end subroutine MPI_SumAll_int_r2

   subroutine MPI_SumAll_int_r3(x)
      implicit none
      !-arguments---------------------------------------
      integer, intent(inout) :: x(:,:,:)
      !-local-variables---------------------------------
      integer :: ierr
      !-------------------------------------------------

#if defined(_HAS_MPI_)
      call MPI_Allreduce(MPI_IN_PLACE, x, size(x), MPI_INTEGER, MPI_SUM, MPI_COMM_WORLD, ierr)
#endif
   end subroutine MPI_SumAll_int_r3

   subroutine MPI_SumAll_real_r0(x)
      implicit none
      !-arguments---------------------------------------
      real(kind=RP), intent(inout) :: x
      !-local-variables---------------------------------
      integer :: ierr
      !-------------------------------------------------

#if defined(_HAS_MPI_)
      call MPI_Allreduce(MPI_IN_PLACE, x, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD, ierr)
#endif
   end subroutine MPI_SumAll_real_r0

   subroutine MPI_SumAll_real_r1(x)
      implicit none
      !-arguments---------------------------------------
      real(kind=RP), intent(inout) :: x(:)
      !-local-variables---------------------------------
      integer :: ierr
      !-------------------------------------------------

#if defined(_HAS_MPI_)
      call MPI_Allreduce(MPI_IN_PLACE, x, size(x), MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD, ierr)
#endif
   end subroutine MPI_SumAll_real_r1

   subroutine MPI_SumAll_real_r2(x)
      implicit none
      !-arguments---------------------------------------
      real(kind=RP), intent(inout) :: x(:,:)
      !-local-variables---------------------------------
      integer :: ierr
      !-------------------------------------------------

#if defined(_HAS_MPI_)
      call MPI_Allreduce(MPI_IN_PLACE, x, size(x), MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD, ierr)
#endif
   end subroutine MPI_SumAll_real_r2

   subroutine MPI_SumAll_real_r3(x)
      implicit none
      !-arguments---------------------------------------
      real(kind=RP), intent(inout) :: x(:,:,:)
      !-local-variables---------------------------------
      integer :: ierr
      !-------------------------------------------------

#if defined(_HAS_MPI_)
      call MPI_Allreduce(MPI_IN_PLACE, x, size(x), MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD, ierr)
#endif
   end subroutine MPI_SumAll_real_r3

   subroutine MPI_MinMax_real_r0(minimum, maximum)
      implicit none
      !-arguments---------------------------------------
      real(kind=RP), intent(inout) :: minimum
      real(kind=RP), intent(inout) :: maximum
      !-local-variables---------------------------------
      integer :: ierr
      integer :: req(2)
      !-------------------------------------------------

#if defined(_HAS_MPI_)
      call MPI_IAllreduce(MPI_IN_PLACE, minimum, 1, MPI_DOUBLE, MPI_MIN, MPI_COMM_WORLD, req(1), ierr)
      call MPI_IAllreduce(MPI_IN_PLACE, maximum, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD, req(2), ierr)
      call MPI_Waitall(2, req, MPI_STATUSES_IGNORE, ierr)
#endif
   end subroutine MPI_MinMax_real_r0

   subroutine MPI_MinMax_real_r1(minimum, maximum)
      implicit none
      !-arguments---------------------------------------
      real(kind=RP), intent(inout) :: minimum(:)
      real(kind=RP), intent(inout) :: maximum(:)
      !-local-variables---------------------------------
      integer :: ierr
      integer :: req(2)
      !-------------------------------------------------

#if defined(_HAS_MPI_)
      call MPI_IAllreduce(MPI_IN_PLACE, minimum, size(minimum), MPI_DOUBLE, MPI_MIN, MPI_COMM_WORLD, req(1), ierr)
      call MPI_IAllreduce(MPI_IN_PLACE, maximum, size(maximum), MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD, req(2), ierr)
      call MPI_Waitall(2, req, MPI_STATUSES_IGNORE, ierr)
#endif
   end subroutine MPI_MinMax_real_r1

   subroutine MPI_MinMax_real_r2(minimum, maximum)
      implicit none
      !-arguments---------------------------------------
      real(kind=RP), intent(inout) :: minimum(:,:)
      real(kind=RP), intent(inout) :: maximum(:,:)
      !-local-variables---------------------------------
      integer :: ierr
      integer :: req(2)
      !-------------------------------------------------

#if defined(_HAS_MPI_)
      call MPI_IAllreduce(MPI_IN_PLACE, minimum, size(minimum), MPI_DOUBLE, MPI_MIN, MPI_COMM_WORLD, req(1), ierr)
      call MPI_IAllreduce(MPI_IN_PLACE, maximum, size(maximum), MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD, req(2), ierr)
      call MPI_Waitall(2, req, MPI_STATUSES_IGNORE, ierr)
#endif
   end subroutine MPI_MinMax_real_r2

   subroutine MPI_MinMax_real_r3(minimum, maximum)
      implicit none
      !-arguments---------------------------------------
      real(kind=RP), intent(inout) :: minimum(:,:,:)
      real(kind=RP), intent(inout) :: maximum(:,:,:)
      !-local-variables---------------------------------
      integer :: ierr
      integer :: req(2)
      !-------------------------------------------------

#if defined(_HAS_MPI_)
      call MPI_IAllreduce(MPI_IN_PLACE, minimum, size(minimum), MPI_DOUBLE, MPI_MIN, MPI_COMM_WORLD, req(1), ierr)
      call MPI_IAllreduce(MPI_IN_PLACE, maximum, size(maximum), MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD, req(2), ierr)
      call MPI_Waitall(2, req, MPI_STATUSES_IGNORE, ierr)
#endif
   end subroutine MPI_MinMax_real_r3

!-----------------------------------------------------------------------------------------------------------------------
! NOTE
!-----------------------------------------------------------------------------------------------------------------------
! This module can be largely simplified with assumed-rank arrays, but they are Fortran 2018 and not all compilers
! work well with this. Take for example the case of the `maximum` variable in `MPI_MinMax`, where the whole construct
! could be simplified to:
!
!   call MPI_IAllreduce(MPI_IN_PLACE, maximum, size(maximum), MPI_DOUBLE, MPI_MIN, MPI_COMM_WORLD, req(1), ierr)
!
! if the argument has the signature:
!
!   real(kind=RP), intent(inout) :: maximum(..)
!-----------------------------------------------------------------------------------------------------------------------

end module MPI_Utilities
