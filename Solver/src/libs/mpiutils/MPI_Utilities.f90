module MPI_Utilities
   use SMConstants
   use MPI_Process_Info
#ifdef _HAS_MPI_
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
      real(kind=RP) :: xloc ! local x
      integer       :: ierr
      !-------------------------------------------------

      if (MPI_Process % doMPIAction) then
         xloc = sum(x**2)
#ifdef _HAS_MPI_
         call mpi_allreduce(xloc, norm, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD, ierr)
#endif
         norm = sqrt(norm)
      else
         norm = sqrt(sum(x*x))
      end if

   end function L2Norm

   function infNorm(x) result(norm)
      implicit none
      !-arguments---------------------------------------
      real(kind=RP), intent(in) :: x(:)
      real(kind=RP)             :: norm
      !-local-variables---------------------------------
      real(kind=RP) :: xloc ! local x
      integer       :: ierr
      !-------------------------------------------------

      norm = maxval(abs(x))

      if (MPI_Process % doMPIAction) then
         xloc = norm
#ifdef _HAS_MPI_
         call mpi_allreduce(xloc, norm, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD, ierr)
#endif
      end if

   end function infNorm

   subroutine MPI_OpAll_int(x, op)
      implicit none
      !-arguments---------------------------------------
      integer, intent(inout) :: x(..)
      integer, intent(in)    :: op
      !-local-variables---------------------------------
      integer :: xloc(size(x)) ! local x
      integer :: ierr
      !-------------------------------------------------

#ifdef _HAS_MPI_
      select rank(x)
      rank(0)
         xloc = x
         call mpi_allreduce(xloc, x, 1, MPI_INTEGER, op, MPI_COMM_WORLD, ierr)
      rank(1)
         xloc = x
         call mpi_allreduce(xloc, x, size(x), MPI_INTEGER, op, MPI_COMM_WORLD, ierr)
      rank(2)
         xloc = pack(x, .true.)
         call mpi_allreduce(xloc, x, size(x), MPI_INTEGER, op, MPI_COMM_WORLD, ierr)
      rank(3)
         xloc = pack(x, .true.)
         call mpi_allreduce(xloc, x, size(x), MPI_INTEGER, op, MPI_COMM_WORLD, ierr)
      end select
#endif
   end subroutine MPI_OpAll_int

   subroutine MPI_OpAll_real(x, op)
      implicit none
      !-arguments---------------------------------------
      real(kind=RP), intent(inout) :: x(..)
      integer,       intent(in)    :: op
      !-local-variables---------------------------------
      real(kind=RP) :: xloc(size(x)) ! local x
      integer       :: ierr
      !-------------------------------------------------

#ifdef _HAS_MPI_
      select rank(x)
      rank(0)
         xloc = x
         call mpi_allreduce(xloc, x, 1, MPI_DOUBLE, op, MPI_COMM_WORLD, ierr)
      rank(1)
         xloc = x
         call mpi_allreduce(xloc, x, size(x), MPI_DOUBLE, op, MPI_COMM_WORLD, ierr)
      rank(2)
         xloc = pack(x, .true.)
         call mpi_allreduce(xloc, x, size(x), MPI_DOUBLE, op, MPI_COMM_WORLD, ierr)
      rank(3)
         xloc = pack(x, .true.)
         call mpi_allreduce(xloc, x, size(x), MPI_DOUBLE, op, MPI_COMM_WORLD, ierr)
      end select
#endif
   end subroutine MPI_OpAll_real

   subroutine MPI_SumAll_int(x)
      implicit none
      !-arguments---------------------------------------
      integer, intent(inout) :: x(..)
      !-------------------------------------------------

#ifdef _HAS_MPI_
      call MPI_OpAll(x, MPI_SUM)
#endif
   end subroutine MPI_SumAll_int

   subroutine MPI_SumAll_real(x)
      implicit none
      !-arguments---------------------------------------
      real(kind=RP), intent(inout) :: x(..)
      !-------------------------------------------------

#ifdef _HAS_MPI_
      call MPI_OpAll(x, MPI_SUM)
#endif
   end subroutine MPI_SumAll_real

   subroutine MPI_MinMax(minimum, maximum)
      implicit none
      !-arguments---------------------------------------
      real(kind=RP), intent(inout) :: minimum(..), maximum(..)
      !-local-variables---------------------------------
      real(kind=RP), allocatable :: lminmax_in(:), lminmax_out(:)
      integer                    :: nmin, nmax, n
      integer                    :: ierr
      !-------------------------------------------------

#ifdef _HAS_MPI_
      select rank(minimum)
      rank(0)
         select rank(maximum)
         rank(0)

            lminmax_in = [-minimum, maximum]
            allocate(lminmax_out(2))

            call mpi_allreduce(lminmax_in, lminmax_out, 2, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD, ierr)

            minimum = -lminmax_out(1)
            maximum = lminmax_out(2)

         end select

      rank(1)
         select rank(maximum)
         rank(1)

            nmin = size(minimum)
            nmax = size(maximum)
            n = nmin + nmax

            lminmax_in = [-minimum, maximum]
            allocate(lminmax_out(n))

            call mpi_allreduce(lminmax_in, lminmax_out, n, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD, ierr)

            minimum = -lminmax_out(1:nmin)
            maximum = lminmax_out(nmin+1:n)

         end select

      end select

#endif
   end subroutine MPI_MinMax

end module MPI_Utilities
