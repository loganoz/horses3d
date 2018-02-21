module GenericMatrixClass
   use SMConstants
   implicit none
   
   private
   public Matrix_t
   
   type Matrix_t
      integer :: NumRows
      contains
         procedure :: construct
         procedure :: Preallocate
         procedure :: Reset
         procedure :: SetColumn
         procedure :: AddToColumn
         procedure :: Shift
         procedure :: ReShift
         procedure :: PreAssembly
         procedure :: Assembly
         procedure :: destruct
   end type Matrix_t
contains
!
!///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
!
   subroutine construct(this,dimPrb,withMPI)
      implicit none
      !---------------------------------------------
      class(Matrix_t) :: this
      integer, intent(in)           :: dimPrb
      logical, optional, intent(in) :: WithMPI
      !---------------------------------------------
      ERROR stop ' :: construct not implemented for current matrix type'
   end subroutine construct
!
!///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
!
   subroutine Preallocate(this, nnz, nnzs)
      implicit none
      !---------------------------------------------
      class(Matrix_t), intent(inout) :: this
      integer, optional, intent(in)  :: nnz
      integer, optional, intent(in)  :: nnzs(:)
      !---------------------------------------------
      ERROR stop ' :: Preallocate not implemented for current matrix type'
   end subroutine Preallocate
!
!///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
!
   subroutine Reset(this)
      implicit none
      !---------------------------------------------
      class(Matrix_t),     intent(inout)     :: this
      !---------------------------------------------
      ERROR stop ' :: Reset not implemented for current matrix type'
   end subroutine Reset
!
!///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
!
   subroutine SetColumn(this,nvalues, irow, icol, values )
      implicit none
      !---------------------------------------------
      class(Matrix_t)            , intent(inout) :: this
      integer                    , intent(in)    :: nvalues
      integer, dimension(:)      , intent(in)    :: irow
      integer                    , intent(in)    :: icol
      real(kind=RP), dimension(:), intent(in)    :: values
      !---------------------------------------------
      ERROR stop ' :: SetColumn not implemented for current matrix type'
   end subroutine SetColumn
!
!///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
! 
   subroutine AddToColumn(this,nvalues, irow, icol, values )
      implicit none
      !---------------------------------------------
      class(Matrix_t), intent(inout)         :: this
      integer, intent(in)                               :: nvalues
      integer, dimension(:), intent(in)                 :: irow
      integer, intent(in)                               :: icol
      real(kind=RP) , dimension(:), intent(in)                 :: values
      !---------------------------------------------
      ERROR stop ' :: AddToColumn not implemented for current matrix type'
   end subroutine AddToColumn
!
!///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
!
   subroutine Shift(this, shiftval)
      implicit none
      !---------------------------------------------
      class(Matrix_t), intent(inout)     :: this
      real(kind=RP),                     intent(in)        :: shiftval
      !---------------------------------------------
      ERROR stop ' :: Shift not implemented for current matrix type'
   end subroutine Shift
!
!///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
!
!  --------------------------------------------------
!  Removes previous shift in order to insert new one 
!              (important when Jacobian is reused)
!  --------------------------------------------------
   subroutine ReShift(this, shiftval)
      implicit none
      !---------------------------------------------
      class(Matrix_t), intent(inout)     :: this
      real(kind=RP),               intent(in)        :: shiftval
      !---------------------------------------------
      ERROR stop ' :: ReShift not implemented for current matrix type'
   end subroutine ReShift
!
!///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
!
   subroutine PreAssembly(this)
      implicit none
      !---------------------------------------------
      class(Matrix_t),     intent(inout)   :: this
      !---------------------------------------------
      ERROR stop ' :: PreAssembly not implemented for current matrix type'
   end subroutine PreAssembly
!
!///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
!
   subroutine Assembly(this,BlockIdx,BlockSize)
      implicit none
      !---------------------------------------------
      class(Matrix_t),     intent(inout)   :: this
      integer, target, optional    ,     intent(in)      :: BlockIdx(:)
      integer, target, optional, intent(in)    :: BlockSize(:)
      !---------------------------------------------
   end subroutine Assembly
!
!///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
!
   subroutine destruct(this)
      implicit none
      !---------------------------------------------
      class(Matrix_t),     intent(inout)     :: this
      !---------------------------------------------
      ERROR stop ' :: destruct not implemented for current matrix type'
   end subroutine destruct
   
end module GenericMatrixClass
