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
         call MPI_Allreduce(MPI_IN_PLACE, x, size(x), MPI_INTEGER, op, MPI_COMM_WORLD, ierr)
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
         call MPI_Allreduce(MPI_IN_PLACE, x, size(x), MPI_DOUBLE, op, MPI_COMM_WORLD, ierr)
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
         call MPI_IAllreduce(MPI_IN_PLACE, minimum, size(minimum), MPI_DOUBLE, MPI_MIN, MPI_COMM_WORLD, req(1), ierr)
         call MPI_IAllreduce(MPI_IN_PLACE, maximum, size(maximum), MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD, req(2), ierr)
         call MPI_Waitall(2, req, MPI_STATUSES_IGNORE, ierr)
      end if
#endif
   end subroutine MPI_MinMax

end module MPI_Utilities