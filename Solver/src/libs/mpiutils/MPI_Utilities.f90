!
!//////////////////////////////////////////////////////
!
!   @File:    MPI_Utilities.f90
!   @Author:  Andrés Rueda (am.rueda@upm.es)
!   @Created: Thu May 16 12:17:25 2019
!   @Last revision date:
!   @Last revision author:
!   @Last revision commit:
!
!//////////////////////////////////////////////////////
!
module MPI_Utilities
   use SMConstants
   use MPI_Process_Info
#ifdef _HAS_MPI_
   use mpi
#endif
   
   private
   public infNorm, L2Norm, MPI_SumAll
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
   
   subroutine MPI_SumAll(x)
      implicit none
      !-arguments---------------------------------------
      real(kind=RP), intent(inout) :: x
      !-local-variables---------------------------------
      real(kind=RP) :: xloc ! local x
      integer       :: ierr
      !-------------------------------------------------
      
#ifdef _HAS_MPI_
      xloc = x
      call mpi_allreduce(xloc, x, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD, ierr)
#endif
   end subroutine MPI_SumAll
   
end module MPI_Utilities
