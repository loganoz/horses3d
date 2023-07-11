#include "Includes.h"
module MPI_Utilities
   use SMConstants
   use MPI_Process_Info
#if defined(_HAS_MPI_)
   use mpi
#endif

   private
   public infNorm, L2Norm, MPI_OpAll, MPI_SumAll, MPI_MinMax

   interface MPI_OpAll
      module procedure MPI_OpAll_int
      module procedure MPI_OpAll_real
   end interface MPI_OpAll

   interface MPI_SumAll
      module procedure MPI_SumAll_int
      module procedure MPI_SumAll_real
   end interface MPI_SumAll
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

   subroutine MPI_OpAll_int(x, op)
      implicit none
      !-arguments---------------------------------------
      integer, intent(inout) :: x(..)
      integer, intent(in)    :: op
      !-local-variables---------------------------------
      integer :: ierr
      !-------------------------------------------------

#if defined(_HAS_MPI_)
      if (MPI_Process % doMPIAction) then
         select rank(x)
         rank(0)
            call MPI_Allreduce(MPI_IN_PLACE, x, 1, MPI_INTEGER, op, MPI_COMM_WORLD, ierr)
         rank(1)
            call MPI_Allreduce(MPI_IN_PLACE, x, size(x), MPI_INTEGER, op, MPI_COMM_WORLD, ierr)
         rank(2)
            call MPI_Allreduce(MPI_IN_PLACE, x, size(x), MPI_INTEGER, op, MPI_COMM_WORLD, ierr)
         rank(3)
            call MPI_Allreduce(MPI_IN_PLACE, x, size(x), MPI_INTEGER, op, MPI_COMM_WORLD, ierr)
         rank(4)
            call MPI_Allreduce(MPI_IN_PLACE, x, size(x), MPI_INTEGER, op, MPI_COMM_WORLD, ierr)
         rank default
            write(STD_OUT,*) 'MPI_OpAll_int only implemented for rank(x) <= 4'
            errorMessage(STD_OUT)
            stop
         end select
      end if
#endif
   end subroutine MPI_OpAll_int

   subroutine MPI_OpAll_real(x, op)
      implicit none
      !-arguments---------------------------------------
      real(kind=RP), intent(inout) :: x(..)
      integer,       intent(in)    :: op
      !-local-variables---------------------------------
      integer       :: ierr
      !-------------------------------------------------

#if defined(_HAS_MPI_)
      if (MPI_Process % doMPIAction) then
         select rank(x)
         rank(0)
            call MPI_Allreduce(MPI_IN_PLACE, x, 1, MPI_DOUBLE, op, MPI_COMM_WORLD, ierr)
         rank(1)
            call MPI_Allreduce(MPI_IN_PLACE, x, size(x), MPI_DOUBLE, op, MPI_COMM_WORLD, ierr)
         rank(2)
            call MPI_Allreduce(MPI_IN_PLACE, x, size(x), MPI_DOUBLE, op, MPI_COMM_WORLD, ierr)
         rank(3)
            call MPI_Allreduce(MPI_IN_PLACE, x, size(x), MPI_DOUBLE, op, MPI_COMM_WORLD, ierr)
         rank(4)
            call MPI_Allreduce(MPI_IN_PLACE, x, size(x), MPI_DOUBLE, op, MPI_COMM_WORLD, ierr)
         rank default
            write(STD_OUT,*) 'MPI_OpAll_real only implemented for rank(x) <= 4'
            errorMessage(STD_OUT)
            stop
         end select
      end if
#endif
   end subroutine MPI_OpAll_real

   subroutine MPI_SumAll_int(x)
      implicit none
      !-arguments---------------------------------------
      integer, intent(inout) :: x(..)
      !-------------------------------------------------

#if defined(_HAS_MPI_)
      call MPI_OpAll(x, MPI_SUM)
#endif
   end subroutine MPI_SumAll_int

   subroutine MPI_SumAll_real(x)
      implicit none
      !-arguments---------------------------------------
      real(kind=RP), intent(inout) :: x(..)
      !-------------------------------------------------

#if defined(_HAS_MPI_)
      call MPI_OpAll(x, MPI_SUM)
#endif
   end subroutine MPI_SumAll_real

   subroutine MPI_MinMax(minimum, maximum)
      implicit none
      !-arguments---------------------------------------
      real(kind=RP), intent(inout) :: minimum(..), maximum(..)
      !-local-variables---------------------------------
      integer :: ierr
      integer :: req(2)
      !-------------------------------------------------

#if defined(_HAS_MPI_)
      if (MPI_Process % doMPIAction) then
         if (rank(minimum) > 4 .or. rank(maximum) > 4) then
            write(STD_OUT,*) 'MPI_MinMax only implemented for rank(x) <= 4'
            errorMessage(STD_OUT)
            stop
         end if

         select rank(minimum)
         rank(0)
            call MPI_IAllreduce(MPI_IN_PLACE, minimum, 1, MPI_DOUBLE, MPI_MIN, MPI_COMM_WORLD, req(1), ierr)
         rank(1)
            call MPI_IAllreduce(MPI_IN_PLACE, minimum, size(minimum), MPI_DOUBLE, MPI_MIN, MPI_COMM_WORLD, req(1), ierr)
         rank(2)
            call MPI_IAllreduce(MPI_IN_PLACE, minimum, size(minimum), MPI_DOUBLE, MPI_MIN, MPI_COMM_WORLD, req(1), ierr)
         rank(3)
            call MPI_IAllreduce(MPI_IN_PLACE, minimum, size(minimum), MPI_DOUBLE, MPI_MIN, MPI_COMM_WORLD, req(1), ierr)
         rank(4)
            call MPI_IAllreduce(MPI_IN_PLACE, minimum, size(minimum), MPI_DOUBLE, MPI_MIN, MPI_COMM_WORLD, req(1), ierr)
         end select

         select rank(maximum)
         rank(0)
            call MPI_IAllreduce(MPI_IN_PLACE, maximum, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD, req(2), ierr)
            call MPI_Waitall(2, req, MPI_STATUSES_IGNORE, ierr)
         rank(1)
            call MPI_IAllreduce(MPI_IN_PLACE, maximum, size(maximum), MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD, req(2), ierr)
            call MPI_Waitall(2, req, MPI_STATUSES_IGNORE, ierr)
         rank(2)
            call MPI_IAllreduce(MPI_IN_PLACE, maximum, size(maximum), MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD, req(2), ierr)
            call MPI_Waitall(2, req, MPI_STATUSES_IGNORE, ierr)
         rank(3)
            call MPI_IAllreduce(MPI_IN_PLACE, maximum, size(maximum), MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD, req(2), ierr)
            call MPI_Waitall(2, req, MPI_STATUSES_IGNORE, ierr)
         rank(4)
            call MPI_IAllreduce(MPI_IN_PLACE, maximum, size(maximum), MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD, req(2), ierr)
            call MPI_Waitall(2, req, MPI_STATUSES_IGNORE, ierr)
         end select
      end if
#endif
   end subroutine MPI_MinMax

!-----------------------------------------------------------------------------------------------------------------------
! NOTE
!-----------------------------------------------------------------------------------------------------------------------
! We have been forced to use `select rank` constructs to avoid compilation issues in several clusters. If the external
! MPI library implements `.mod` files that support assumed-rank variables, the code can be simplified. Take for example
! the case of the `maximum` variable in `MPI_MinMax`, where the whole construct:
!
!   select rank(maximum)
!   rank(0)
!     call MPI_IAllreduce(MPI_IN_PLACE, maximum, 1, MPI_DOUBLE, MPI_MIN, MPI_COMM_WORLD, req(1), ierr)
!   rank(1)
!     call MPI_IAllreduce(MPI_IN_PLACE, maximum, size(maximum), MPI_DOUBLE, MPI_MIN, MPI_COMM_WORLD, req(1), ierr)
!   ...
!   end select
!
! can be simplified to:
!
!   call MPI_IAllreduce(MPI_IN_PLACE, maximum, size(maximum), MPI_DOUBLE, MPI_MIN, MPI_COMM_WORLD, req(1), ierr)
!-----------------------------------------------------------------------------------------------------------------------

end module MPI_Utilities
